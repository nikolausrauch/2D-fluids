#include "iisph.h"

#include "utility/perfmonitor.h"
#include "utility/particlehelper.h"

#include <cassert>
#include <algorithm>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*128)

Particle::Particle(const glm::vec2 &pos)
    : position(pos), velocity{{0.0f, 0.0f}, {0.0f, 0.0f}},
      force(0.0f),
      density(0.0f), densityAdv(0.0f), pressure(0.0f),
      d_ii(0.0f, 0.0f), sum_d_ij_p_j(0.0f, 0.0f), a_ii(0.0f), pressure_i{0.0f, 0.0f}
{

}

GhostParticle::GhostParticle(const glm::vec2& pos)
    : position(pos)
{

}

IISPH::IISPH()
    : timeStep(0.0064f), stepsPerFrame(4),
      mass(0.006f), radiusParticle(0.05), radiusKernel(kernel::radiusMultiplier*radiusParticle),
      restDensity(1.0f), viscocityConstant(0.1f),
      numIterations(7),
      kernelDensity(radiusKernel), kernelPressure(radiusKernel), kernelViscocity(radiusKernel),
      gravity(0.0f, -9.81f),
      fluidNNsearch(PARTICLE_INIT), boundaryNNsearch(PARTICLE_INIT)
{
    fluidParticles.reserve(PARTICLE_INIT);
    autotuneParameter();
}

void IISPH::create(const Scene& desc)
{
    clear();

    for(const auto& box : desc.boxes())
    {
        if(box.type == Scene::eType::BOUNDARY)
        {
            auto min = box.min + 0.5f*glm::vec2{radiusParticle, radiusParticle};
            auto max = box.max - 0.5f*glm::vec2{radiusParticle, radiusParticle};
            auto sampled_pos = helper::randomPositions(radiusParticle, min, max);
            boundaryParticles.insert(boundaryParticles.end(), sampled_pos.begin(), sampled_pos.end());
        }
        else if(box.type == Scene::eType::FLUID_BODY)
        {
            auto min = box.min + glm::vec2{radiusParticle, radiusParticle};
            auto max = box.max - glm::vec2{radiusParticle, radiusParticle};
            auto sampled_pos = helper::randomPositions(2.0f*radiusParticle, min, max);
            fluidParticles.insert(fluidParticles.end(), sampled_pos.begin(), sampled_pos.end());
        }
    }

    for(const auto& circle : desc.circles())
    {
        if(circle.type == Scene::eType::BOUNDARY)
        {
            auto sampled_pos = helper::randomPositions(radiusParticle, circle.center, circle.radius);
            boundaryParticles.insert(boundaryParticles.end(), sampled_pos.begin(), sampled_pos.end());
        }
        else if(circle.type == Scene::eType::FLUID_BODY)
        {
            auto sampled_pos = helper::randomPositions(2.0f*radiusParticle, circle.center, circle.radius);
            fluidParticles.insert(fluidParticles.end(), sampled_pos.begin(), sampled_pos.end());
        }
    }

    /* assume stationary boundary particles */
    boundaryNNsearch.fillGrid(boundaryParticles, radiusKernel);
}

Particle& IISPH::createParticle(const glm::vec2 &pos)
{
    assert(particles.size() < PARTICLE_INIT);

    fluidParticles.emplace_back(pos);
    return fluidParticles.back();
}

void IISPH::createParticles(const std::vector<glm::vec2>& pos)
{
    fluidParticles.insert(fluidParticles.end(), pos.begin(), pos.end());
}

GhostParticle& IISPH::createGhostParticle(const glm::vec2& pos)
{
    assert(particles.size() < PARTICLE_INIT);

    boundaryParticles.emplace_back(pos);
    return boundaryParticles.back();
}

unsigned int IISPH::index(const Particle& p) const
{
    return std::distance(fluidParticles.data(), &p);
}

unsigned int IISPH::index(const GhostParticle& p) const
{
    return std::distance(boundaryParticles.data(), &p);
}

void IISPH::stepFillNNsearch()
{
    fluidNNsearch.clear();
    fluidNNsearch.fillGrid(fluidParticles, radiusKernel);
}

void IISPH::stepPredictVelocity(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        fluidNNsearch.updateNeighbour(index(p_i), fluidParticles, radiusKernel);
        boundaryNNsearch.updateGhostNeighbour(index(p_i), p_i.position, boundaryParticles, radiusKernel);

        /* density */
        p_i.density = mass * kernelDensity(0, 0);

        /* fluid particles */
        {
            const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto d = glm::length(p_i.position - fluidParticles[j].position);
                p_i.density += mass * kernelDensity(d, d*d);
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto d = glm::length(p_i.position - boundaryParticles[j].position);
                p_i.density += mass * kernelDensity(d, d*d);
            });
        }


        /* gravity */
        p_i.force = mass * gravity;

        /* viscocity */
        const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
        {
            const auto& p_j = fluidParticles[j];

            p_i.force += mass * mass * viscocityConstant * radiusKernel * 88.5f
                    * glm::dot(p_i.velocity[0] - p_j.velocity[0], p_i.position - p_j.position)
                    / (glm::dot(p_i.position - p_j.position, p_i.position - p_j.position) + 0.1f * radiusKernel * radiusKernel)
                    * kernelViscocity.gradient(p_i.position, p_j.position);
        });


        /* predict velocity */
        p_i.velocity[1] = p_i.velocity[0] + dt / mass * p_i.force;


        /* precompute some of the linear system constants */
        p_i.d_ii = {0.0f, 0.0f};

        /* fluid particles */
        {
            const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                p_i.d_ii -= dt*dt* mass / (p_i.density*p_i.density) * kernelPressure.gradient(p_i.position, p_j.position);
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                p_i.d_ii -= dt*dt* mass / (p_i.density*p_i.density) * kernelPressure.gradient(p_i.position, p_j.position);
            });
        }
    }
}

void IISPH::stepAdvectDensity(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        profile_sample(advection);

        /* initialize unknown pressure from prev. solution */
        p_i.pressure_i[0] = 0.5f * p_i.pressure;

        /* density */
        p_i.densityAdv = p_i.density;

        /* fluid particles */
        {
            const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                p_i.densityAdv += dt * mass * glm::dot(p_i.velocity[1] - p_j.velocity[1], kernelDensity.gradient(p_i.position, p_j.position));
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                /* assumes stationary boundaries -> otherwise p_i.velocity[1] - p_j.velocity */
                const auto& p_j = boundaryParticles[j];
                p_i.densityAdv += dt * mass * glm::dot(p_i.velocity[1], kernelDensity.gradient(p_i.position, p_j.position));
            });
        }


        /* compute  diagonal entries a_ii */
        p_i.a_ii = 0.0f;

        /* fluid particles */
        {
            const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                auto d_ii_d_ji = (p_i.d_ii - (-dt*dt* mass / (p_i.density*p_i.density) * kernelPressure.gradient(p_j.position, p_i.position)) );
                p_i.a_ii += mass * glm::dot(d_ii_d_ji, kernelPressure.gradient(p_i.position, p_j.position));
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                auto d_ii_d_ji = (p_i.d_ii - (-dt*dt* mass / (p_i.density*p_i.density) * kernelPressure.gradient(p_j.position, p_i.position)) );
                p_i.a_ii += mass * glm::dot(d_ii_d_ji, kernelPressure.gradient(p_i.position, p_j.position));
            });
        }
    }
}

void IISPH::stepComputeCoeff(float dt, int idx, int idx_2)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.sum_d_ij_p_j = {0.0f, 0.0f};

        const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
        {
            const auto& p_j = fluidParticles[j];
            p_i.sum_d_ij_p_j -= dt*dt * mass / (p_j.density*p_j.density) * p_j.pressure_i[idx] * kernelPressure.gradient(p_i.position, p_j.position);
        });
    }
}

void IISPH::stepUpdatePressure(float dt, int idx, int idx_2)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        float sum_j = 0.0f;
        float sum_b = 0.0f;

        /* fluid particles */
        {
            const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                auto sum_d_jk = p_j.sum_d_ij_p_j - (-dt*dt* mass / (p_i.density*p_i.density) * p_i.pressure_i[idx] * kernelPressure.gradient(p_j.position, p_i.position));
                sum_j += mass * glm::dot(p_i.sum_d_ij_p_j - p_j.d_ii * p_j.pressure_i[idx] - sum_d_jk, kernelPressure.gradient(p_i.position, p_j.position));
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                sum_b += mass * glm::dot(p_i.sum_d_ij_p_j, kernelPressure.gradient(p_i.position, p_j.position));
            });
        }

        float relax = 0.5f;
        if(std::fabs(p_i.a_ii) > 1e-9f)
        {
            p_i.pressure_i[idx_2] = (1.0f - relax) * p_i.pressure_i[idx] + relax / p_i.a_ii * (restDensity - p_i.densityAdv - sum_j - sum_b);
        }
        else
        {
            p_i.pressure_i[idx_2] = 0.0f;
        }

        p_i.pressure_i[idx_2] = std::max(p_i.pressure_i[idx_2], 0.0f);
        p_i.pressure = p_i.pressure_i[idx_2];
    }
}

void IISPH::stepPressureForce()
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {

        p_i.force = {0.0f, 0.0f};

        /* fluid particles */
        {
            const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];

                /* pressure gradient */
                p_i.force -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_j.pressure / (p_j.density*p_j.density))
                        * kernelPressure.gradient(p_i.position, p_j.position);
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                p_i.force -= mass * mass * (p_i.pressure / (p_i.density*p_i.density))
                        * kernelPressure.gradient(p_i.position, p_j.position);
            });
        }
    }
}

void IISPH::stepIntegrate(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.velocity[0] = p_i.velocity[1] + dt * p_i.force / mass;
        p_i.position += dt * p_i.velocity[0];
    };
}

void IISPH::update()
{
    PerfMonitor::instance().start_frame();
    profile_sample(iisph_total);

    for(int i = 0; i < stepsPerFrame; i++)
    {
        update(timeStep);
    }
}

void IISPH::update(float dt)
{
    {
        profile_sample(iisph_nnfill);
        stepFillNNsearch();
    }

    {
        profile_sample(iisph_predvel);
        stepPredictVelocity(dt);
    }

    {
        profile_sample(iisph_advectdensity);
        stepAdvectDensity(dt);
    }

    for(int i = 0; i < numIterations; i++)
    {
        profile_sample(iisph_pressuresolve);

        int p_idx = i % 2;
        int p_idx_1 = (i+1) % 2;

        {
            profile_sample(iisph_coeff);
            stepComputeCoeff(dt, p_idx, p_idx_1);
        }

        {
            profile_sample(iisph_updatepress);
            stepUpdatePressure(dt, p_idx, p_idx_1);
        }
    }

    {
        profile_sample(iisph_pressureforce);
        stepPressureForce();
    }

    {
        profile_sample(iisph_integrate);
        stepIntegrate(dt);
    }
}

void IISPH::updateRadius(float radius)
{
    radiusParticle = radius;
    radiusKernel = kernel::radiusMultiplier * radius;

    kernelDensity = kernel::std(radiusKernel);
    kernelPressure = kernel::spiky(radiusKernel);
    kernelViscocity = kernel::std(radiusKernel);

    autotuneParameter();
}

void IISPH::autotuneParameter()
{
    auto sampledPositions = helper::randomPositions(2.0f*radiusParticle, {-2.0f*radiusKernel, -2.0f*radiusKernel}, {2.0f*radiusKernel, 2.0f*radiusKernel});
    float densityEstimate = 0.0f;

    for(const auto& p : sampledPositions)
    {
        float d = glm::length(p);
        densityEstimate += kernelDensity(d, d*d);
    }

    mass = restDensity / densityEstimate;
}

void IISPH::clear()
{
    fluidParticles.clear();
    boundaryParticles.clear();

    fluidNNsearch.clear();
    boundaryNNsearch.clear();
}
