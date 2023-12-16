#include "pcisph.h"

#include <utility/perfmonitor.h>

#include <cassert>
#include <algorithm>
#include <iostream>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*256)
#define RAND_0_1 (static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX))


Particle::Particle(const glm::vec2 &pos)
    : position{pos, {0.0f,0.0f}}, velocity{{0.0, 0.0}, {0.0f, 0.0f}}, forceNoPress(0.0, 0.0), forcePress(0.0f, 0.0f)
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


PCISPH::PCISPH()
    : timeStep(0.0032), numPerFrame(8),
      mass(0.006), radiusParticle(0.05), radiusKernel(4.0f*radiusParticle),
      restDensity(1.0), viscocityConstant(0.1f),
      numIterations(3), errorThresh(0.01f), stiffnessGradFact(1.0f),
      gravity(0, -9.81),
      densityKernel(radiusKernel), pressureKernel(radiusKernel), viscosityKernel(radiusKernel),
      nnSearch(PARTICLE_INIT), nnSearchBoundary(PARTICLE_INIT)
{
    particles.reserve(PARTICLE_INIT);
    autotuneParams();
}

Particle &PCISPH::createParticle(const glm::vec2 &pos)
{
    assert(particles.size() < PARTICLE_INIT);

    auto& p = particles.emplace_back(pos);
    p.density = restDensity;
    return particles.back();
}

GhostParticle& PCISPH::createGhostParticle(const glm::vec2& pos)
{
    assert(particles.size() < PARTICLE_INIT);

    particlesBoundary.emplace_back(pos);
    return particlesBoundary.back();
}

Particle &PCISPH::particle(unsigned int index)
{
    assert(index < particles.size());

    return particles[index];
}

const Particle &PCISPH::particle(unsigned int index) const
{
    assert(index < particles.size());

    return particles[index];
}

unsigned int PCISPH::index(const Particle &p) const
{
    return std::distance(particles.data(), &p);
}

unsigned int PCISPH::index(const Particle *p) const
{
    return std::distance(particles.data(), p);
}

void PCISPH::boundary(const glm::vec2 &pos, const glm::vec2 &normal, const glm::vec2& min, const glm::vec2& max)
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

void PCISPH::clear()
{
    particles.clear();
    particlesBoundary.clear();

    nnSearch.clear();
    nnSearchBoundary.clear();
}

void PCISPH::update()
{
    PerfMonitor::instance().start_frame();
    profile_sample(total);

    for(int i = 0; i < numPerFrame; i++)
    {
        update(timeStep);
    }
}

void PCISPH::update(float dt)
{
    {
        profile_sample(nnfill);
        nnSearch.clear();
        nnSearch.fillGrid(particles, radiusKernel);
    }

    {
        profile_sample(extforce);

        #pragma omp parallel for schedule(static)
        for(auto& p_i : particles)
        {
            nnSearch.updateNeighbour(std::distance(particles.data(), &p_i), particles, radiusKernel);
            nnSearchBoundary.updateGhostNeighbour(std::distance(particles.data(), &p_i), p_i.position[0], particlesBoundary, radiusKernel);

            /* gravity */
            p_i.forceNoPress = mass * gravity;

            /* viscocity */
            const auto& neighbours = nnSearch.neighbors(p_i, particles);
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = particles[j];

                /* viscocity */
                p_i.forceNoPress += mass * mass * viscocityConstant * radiusKernel * 88.5f
                        * glm::dot(p_i.velocity[0] - p_j.velocity[0], p_i.position[0] - p_j.position[0])
                        / (glm::dot(p_i.position[0] - p_j.position[0], p_i.position[0] - p_j.position[0]) + 0.1f * radiusKernel * radiusKernel)
                        * viscosityKernel.gradient(p_i.position[0], p_j.position[0]);

            });

            /* prepare for prediction-correction step */
            p_i.forcePress = {0.0f, 0.0f};
            p_i.pressure = 0.0f; 
        }

        {
            profile_sample(predcorr);

            /* correction scaling factor */
            float beta = 2.0f * (mass*mass * dt*dt / (restDensity*restDensity));
            float k = -1.0f / ( beta * -stiffnessGradFact );

            float avg_error = 1.0f;
            for(int i = 0; i < numIterations && avg_error > errorThresh; i++)
            {
                /* predict velocity and position */
                #pragma omp parallel for schedule(static)
                for(auto& p_i : particles)
                {
                    p_i.velocity[1] = p_i.velocity[0] + dt / mass * (p_i.forceNoPress + p_i.forcePress);
                    p_i.position[1] = p_i.position[0] + dt * p_i.velocity[1];
                }

                /* compute density at predicted positions */
                #pragma omp parallel for schedule(static)
                for(auto& p : particles)
                {
                    p.density = mass * densityKernel(0, 0);

                    /* fluid particles */
                    {
                        const auto& neighbours = nnSearch.neighbors(p, particles);
                        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                        {
                            const auto d = glm::length(p.position[1] - particles[j].position[1]);
                            p.density += mass * densityKernel(d, d*d);
                        });
                    }

                    /* boundary particles */
                    {
                        const auto& neighbours = nnSearchBoundary.neighbors(std::distance(particles.data(), &p));
                        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                        {
                            const auto d = glm::length(p.position[1] - particlesBoundary[j].position);
                            p.density += mass * restDensity * densityKernel(d, d*d);
                        });
                    }

                    p.densityError = p.density - restDensity;
                    p.pressure += k * p.densityError;
                    p.pressure = std::max(p.pressure, 0.0f);
                }

                /* compute avg error */
                avg_error = 0.0f;
                #pragma omp parallel for schedule(static) reduction(+:avg_error)
                for(auto& p : particles)
                {
                    avg_error += std::max(p.densityError, 0.0f);
                }
                avg_error = avg_error / static_cast<float>(particles.size());

                /* compute pressure gradient */
                #pragma omp parallel for schedule(static)
                for(auto& p_i : particles)
                {
                    p_i.forcePress = {0.0f, 0.0f};

                    /* fluid particles */
                    {
                        const auto& neighbours = nnSearch.neighbors(p_i, particles);
                        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                        {
                            const auto& p_j = particles[j];
                            p_i.forcePress -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_j.pressure / (p_j.density*p_j.density))
                            * pressureKernel.gradient(p_i.position[0], p_j.position[0]);
                        });
                    }

                    /* boundary particles */
                    {
                        const auto& neighbours = nnSearchBoundary.neighbors(std::distance(particles.data(), &p_i));
                        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                        {
                            const auto& p_j = particlesBoundary[j];
                            p_i.forcePress -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_i.pressure / (restDensity*restDensity))
                                    * pressureKernel.gradient(p_i.position[0], p_j.position);
                        });
                    }
                }
            }
        }

        {
            profile_sample(integrate);

            /* final integration */
            #pragma omp parallel for schedule(static)
            for(auto& p_i : particles)
            {
                p_i.velocity[0] = p_i.velocity[0] + dt / mass * (p_i.forceNoPress + p_i.forcePress);
                p_i.position[0] = p_i.position[0] + dt * p_i.velocity[0];
            }
        }

    }
}

double PCISPH::particleMass() const
{
    return mass;
}

void PCISPH::particleMass(float m)
{
    mass = m;
}

float PCISPH::particleRadius() const
{
    return radiusParticle;
}

void PCISPH::particleRadius(float r)
{
    radiusParticle = r;
    kernelRadius(4.0f*radiusParticle);
}

void PCISPH::kernelRadius(float h)
{
    radiusKernel = h;

    densityKernel = kernel::std(h);
    pressureKernel = kernel::spiky(h);

    autotuneParams();
}

float PCISPH::kernelRadius()
{
    return radiusKernel;
}

void PCISPH::autotuneParams()
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

void PCISPH::computeInitialDensity()
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
