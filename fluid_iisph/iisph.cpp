#include "iisph.h"

#include <utility/perfmonitor.h>

#include <cassert>
#include <algorithm>
#include <iostream>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*256)
#define RAND_0_1 (static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX))


Particle::Particle(const glm::vec2 &pos)
    : position{pos, {0.0f,0.0f}}, velocity{{0.0, 0.0}, {0.0f, 0.0f}}, forceNoPress(0.0, 0.0), forcePress(0.0f, 0.0f),
      density(0.0f), densityAdv(0.0f), pressure(0.0f),
      d_ii(0.0f, 0.0f), sum_d_ij_p_j(0.0f, 0.0f), a_ii(0.0f), p_i{0.0f, 0.0f}
{

}

GhostParticle::GhostParticle(const glm::vec2& pos)
    : position(pos)
{

}

Boundary::Boundary(const glm::vec2 &pos, const glm::vec2 &norm, const glm::vec2& min, const glm::vec2& max)
    : position(pos), normal(norm), min(min), max(max)
{

}

bool Boundary::inside(const glm::vec2 &pos)
{
    return penetration(pos) > 0;
}

float Boundary::penetration(const glm::vec2 &pos)
{
    float depth = glm::dot( (position - pos), normal);
    return depth;
}


IISPH::IISPH()
    : timeStep(0.0064), numPerFrame(4),
      mass(0.006), radiusParticle(0.05), radiusKernel(4.0f*radiusParticle),
      restDensity(1.0), viscocityConstant(0.1f),
      numIterations(7), errorThresh(0.01f), stiffnessGradFact(1.0f),
      gravity(0, -9.81),
      densityKernel(radiusKernel), pressureKernel(radiusKernel), viscosityKernel(radiusKernel),
      nnSearch(PARTICLE_INIT), nnSearchBoundary(PARTICLE_INIT)
{
    particles.reserve(PARTICLE_INIT);
    autotuneParams();
}

Particle &IISPH::createParticle(const glm::vec2 &pos)
{
    assert(particles.size() < PARTICLE_INIT);

    auto& p = particles.emplace_back(pos);
    p.density = restDensity;
    return particles.back();
}

GhostParticle& IISPH::createGhostParticle(const glm::vec2& pos)
{
    assert(particles.size() < PARTICLE_INIT);

    particlesBoundary.emplace_back(pos);
    return particlesBoundary.back();
}

Particle &IISPH::particle(unsigned int index)
{
    assert(index < particles.size());

    return particles[index];
}

const Particle &IISPH::particle(unsigned int index) const
{
    assert(index < particles.size());

    return particles[index];
}

unsigned int IISPH::index(const Particle &p) const
{
    return std::distance(particles.data(), &p);
}

unsigned int IISPH::index(const Particle *p) const
{
    return std::distance(particles.data(), p);
}

void IISPH::boundary(const glm::vec2 &pos, const glm::vec2 &normal, const glm::vec2& min, const glm::vec2& max)
{
    boundaries.emplace_back(pos, normal, min, max);

    auto samples = helper::randomPositions(radiusParticle, min, max);
    for(const auto& p : samples)
    {
        createGhostParticle(p);
    }

    nnSearchBoundary.clear();
    nnSearchBoundary.fillGrid(particlesBoundary, radiusKernel);
}

void IISPH::clear()
{
    particles.clear();
    particlesBoundary.clear();

    nnSearch.clear();
    nnSearchBoundary.clear();
}

void IISPH::update()
{
    PerfMonitor::instance().start_frame();
    profile_sample(total);

    for(int i = 0; i < numPerFrame; i++)
    {
        update(timeStep);
    }
}

void IISPH::update(float dt)
{
    {
        profile_sample(nnfill);
        nnSearch.clear();
        nnSearch.fillGrid(particles, radiusKernel);
    }

    {
        profile_sample(predvel);

        #pragma omp parallel for schedule(static)
        for(auto& p_i : particles)
        {
            nnSearch.updateNeighbour(index(p_i), particles, radiusKernel);
            nnSearchBoundary.updateGhostNeighbour(index(p_i), p_i.position[0], particlesBoundary, radiusKernel);

            /* density */
            p_i.density = mass * densityKernel(0, 0);

            /* fluid particles */
            {
                const auto& neighbours = nnSearch.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto d = glm::length(p_i.position[0] - particles[j].position[0]);
                    p_i.density += mass * densityKernel(d, d*d);
                });
            }

            /* boundary particles */
            {
                const auto& neighbours = nnSearchBoundary.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto d = glm::length(p_i.position[0] - particlesBoundary[j].position);
                    p_i.density += mass * densityKernel(d, d*d);
                });
            }


            /* gravity */
            p_i.forceNoPress = mass * gravity;

            /* viscocity */
            const auto& neighbours = nnSearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = particles[j];

                p_i.forceNoPress += mass * mass * viscocityConstant * radiusKernel * 88.5f
                        * glm::dot(p_i.velocity[0] - p_j.velocity[0], p_i.position[0] - p_j.position[0])
                        / (glm::dot(p_i.position[0] - p_j.position[0], p_i.position[0] - p_j.position[0]) + 0.1f * radiusKernel * radiusKernel)
                        * viscosityKernel.gradient(p_i.position[0], p_j.position[0]);
            });


            /* predict velocity */
            p_i.velocity[1] = p_i.velocity[0] + dt / mass * p_i.forceNoPress;


            /* precompute some of the linear system constants */
            p_i.d_ii = {0.0f, 0.0f};
            /* fluid particles */
            {
                const auto& neighbours = nnSearch.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particles[j];
                    p_i.d_ii -= dt*dt* mass / (p_i.density*p_i.density) * pressureKernel.gradient(p_i.position[0], p_j.position[0]);
                });
            }

            /* boundary particles */
            {
                const auto& neighbours = nnSearchBoundary.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particlesBoundary[j];
                    p_i.d_ii -= dt*dt* mass / (p_i.density*p_i.density) * pressureKernel.gradient(p_i.position[0], p_j.position);
                });
            }
        }

        /* advected density and precompute diagonal entries of linear system */
        #pragma omp parallel for schedule(static)
        for(auto& p_i : particles)
        {
            profile_sample(advection);

            /* initialize unknown pressure from prev. solution */
            p_i.p_i[0] = 0.5f * p_i.pressure;

            /* density */
            p_i.densityAdv = p_i.density;

            /* fluid particles */
            {
                const auto& neighbours = nnSearch.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particles[j];
                    p_i.densityAdv += dt * mass * glm::dot(p_i.velocity[1] - p_j.velocity[1], densityKernel.gradient(p_i.position[0], p_j.position[0]));
                });
            }

            /* boundary particles */
            {
                const auto& neighbours = nnSearchBoundary.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    /* assumes stationary boundaries -> otherwise p_i.velocity[1] - p_j.velocity */
                    const auto& p_j = particlesBoundary[j];
                    p_i.densityAdv += dt * mass * glm::dot(p_i.velocity[1], densityKernel.gradient(p_i.position[0], p_j.position));
                });
            }


            /* compute  diagonal entries a_ii */
            p_i.a_ii = 0.0f;

            /* fluid particles */
            {
                const auto& neighbours = nnSearch.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particles[j];
                    auto d_ii_d_ji = (p_i.d_ii - (-dt*dt* mass / (p_i.density*p_i.density) * pressureKernel.gradient(p_j.position[0], p_i.position[0])) );
                    p_i.a_ii += mass * glm::dot(d_ii_d_ji, pressureKernel.gradient(p_i.position[0], p_j.position[0]));
                });
            }

            /* boundary particles */
            {
                const auto& neighbours = nnSearchBoundary.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particlesBoundary[j];
                    auto d_ii_d_ji = (p_i.d_ii - (-dt*dt* mass / (p_i.density*p_i.density) * pressureKernel.gradient(p_j.position, p_i.position[0])) );
                    p_i.a_ii += mass * glm::dot(d_ii_d_ji, pressureKernel.gradient(p_i.position[0], p_j.position));
                });
            }
        }
    }

    {
        profile_sample(pressuresolve);

        /* pressure solver with relaxed Jacobi */
        for(int i = 0; i < numIterations; i++)
        {
            int p_idx = i % 2;
            int p_idx_1 = (i+1) % 2;

            /* compute coefficients */
            #pragma omp parallel for schedule(static)
            for(auto& p_i : particles)
            {
                p_i.sum_d_ij_p_j = {0.0f, 0.0f};

                /* fluid particles */
                {
                    const auto& neighbours = nnSearch.neighbors(index(p_i));
                    std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                    {
                        const auto& p_j = particles[j];
                        p_i.sum_d_ij_p_j -= dt*dt * mass / (p_j.density*p_j.density) * p_j.p_i[p_idx] * pressureKernel.gradient(p_i.position[0], p_j.position[0]);
                    });
                }
            }

            /* compute updated pressure */
            #pragma omp parallel for schedule(static)
            for(auto& p_i : particles)
            {
                float sum_j = 0.0f;
                float sum_b = 0.0f;

                /* fluid particles */
                {
                    const auto& neighbours = nnSearch.neighbors(index(p_i));
                    std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                    {
                        const auto& p_j = particles[j];
                        auto sum_d_jk = p_j.sum_d_ij_p_j - (-dt*dt* mass / (p_i.density*p_i.density) * p_i.p_i[p_idx] * pressureKernel.gradient(p_j.position[0], p_i.position[0]));
                        sum_j += mass * glm::dot(p_i.sum_d_ij_p_j - p_j.d_ii * p_j.p_i[p_idx] - sum_d_jk, pressureKernel.gradient(p_i.position[0], p_j.position[0]));
                    });
                }

                /* boundary particles */
                {
                    const auto& neighbours = nnSearchBoundary.neighbors(index(p_i));
                    std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                    {
                        const auto& p_j = particlesBoundary[j];
                        sum_b += mass * glm::dot(p_i.sum_d_ij_p_j, pressureKernel.gradient(p_i.position[0], p_j.position));
                    });
                }

                float relax = 0.5f;

                if(std::fabs(p_i.a_ii) > 1e-9f)
                {
                    p_i.p_i[p_idx_1] = (1.0f - relax) * p_i.p_i[p_idx] + relax / p_i.a_ii * (restDensity - p_i.densityAdv - sum_j - sum_b);
                }
                else
                {
                    p_i.p_i[p_idx_1] = 0.0f;
                }

                p_i.p_i[p_idx_1] = std::max(p_i.p_i[p_idx_1], 0.0f);
                p_i.pressure = p_i.p_i[p_idx_1];
            }
        }
    }


    {
        profile_sample(integrate);

        #pragma omp parallel for schedule(static)
        for(auto& p_i : particles)
        {
            p_i.forcePress = {0.0f, 0.0f};

            /* fluid particles */
            {
                const auto& neighbours = nnSearch.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particles[j];

                    /* pressure gradient */
                    p_i.forcePress -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_j.pressure / (p_j.density*p_j.density))
                            * pressureKernel.gradient(p_i.position[0], p_j.position[0]);
                });
            }

            /* boundary particles */
            {
                const auto& neighbours = nnSearchBoundary.neighbors(index(p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particlesBoundary[j];
                    p_i.forcePress -= mass * mass * (p_i.pressure / (p_i.density*p_i.density))
                            * pressureKernel.gradient(p_i.position[0], p_j.position);
                });
            }
        };

        #pragma omp parallel for schedule(static)
        for(auto& p_i : particles)
        {
            p_i.velocity[0] = p_i.velocity[1] + dt * p_i.forcePress / mass;
            p_i.position[0] += dt * p_i.velocity[0];
        };
    }
}

double IISPH::particleMass() const
{
    return mass;
}

void IISPH::particleMass(float m)
{
    mass = m;
}

float IISPH::particleRadius() const
{
    return radiusParticle;
}

void IISPH::particleRadius(float r)
{
    radiusParticle = r;
    kernelRadius(4.0f*radiusParticle);
}

void IISPH::kernelRadius(float h)
{
    radiusKernel = h;

    densityKernel = kernel::std(h);
    pressureKernel = kernel::spiky(h);

    autotuneParams();
}

float IISPH::kernelRadius()
{
    return radiusKernel;
}

void IISPH::autotuneParams()
{
    float fact = 2.0f;
    auto sampledPositions = helper::randomPositions(2.0f*radiusParticle, {-fact*radiusKernel, -fact*radiusKernel}, {fact*radiusKernel, fact*radiusKernel});

    unsigned int count = 0;
    float densityEstimate = 0.0f;
    glm::vec2 p_grad = {0.0f, 0.0f};
    glm::vec2 d_grad = {0.0f, 0.0f};
    float grad_dot_sum = 0.0f;

    for(const auto& p : sampledPositions)
    {
        float d = glm::length(p);
        densityEstimate += densityKernel(d, d*d);

        if(densityKernel(d, d*d) > 0) { count++; }

        if(glm::length2(p) < radiusParticle*radiusParticle*0.98) continue;
        if(glm::length(p) >= radiusKernel) continue;

        auto d_wij = pressureKernel.gradient({0.0f, 0.0f}, p);
        auto p_wij = pressureKernel.gradient({0.0f, 0.0f}, p);
        p_grad += p_wij;
        d_grad += d_wij;
        grad_dot_sum += glm::dot(p_wij, d_wij);
    }

    mass = restDensity / densityEstimate;
    stiffnessGradFact = glm::dot(p_grad, d_grad) + grad_dot_sum;

    std::cout << "Particle in range: " << count << std::endl;
    std::cout << "Density Estimate: " << densityEstimate * mass << std::endl;
    std::cout << "Autotune Mass: " << mass << std::endl;
    std::cout << "Gradient Stiffness Factor: " << stiffnessGradFact << std::endl;
}

void IISPH::computeInitialDensity()
{
    nnSearch.clear();
    nnSearch.fillGrid(particles, radiusKernel);

    #pragma omp parallel for schedule(static)
    for(auto& p : particles)
    {
        nnSearch.updateNeighbour(std::distance(particles.data(), &p), particles, radiusKernel);
        nnSearchBoundary.updateGhostNeighbour(std::distance(particles.data(), &p), p.position[0], particlesBoundary, radiusKernel);

        p.density = 0.0f;

        /* fluid particles */
        {
            const auto& neighbours = nnSearch.neighbors(p, particles);
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto d = glm::length(p.position[0] - particles[j].position[0]);
                p.density += mass * densityKernel(d, d*d);
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = nnSearchBoundary.neighbors(std::distance(particles.data(), &p));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto d = glm::length(p.position[0] - particlesBoundary[j].position);
                p.density += mass * restDensity * densityKernel(d, d*d);
            });
        }
    }
}

std::vector<glm::vec2> helper::randomPositions(float spacing, const glm::vec2 &min, const glm::vec2 &max, float jitterNoise)
{
    std::srand(std::time(nullptr));

    std::vector<glm::vec2> positions;
    positions.reserve(1024);

    for(float x = min.x; x <= max.x; x += spacing)
    {
        for(float y = min.y; y <= max.y; y += spacing)
        {

            glm::vec2 pos = {x, y};

            pos.x += ((RAND_0_1 <= 0.5) ? -1.0 : 1.0) * jitterNoise * spacing * RAND_0_1;
            pos.y += ((RAND_0_1 <= 0.5) ? -1.0 : 1.0) * jitterNoise * spacing * RAND_0_1;

            positions.emplace_back(pos);
        }
    }

    return positions;
}
