#include "pbfsystem.h"

#include <utility/perfmonitor.h>

#include <cmath>
#include <ctime>
#include <cassert>
#include <algorithm>
#include <iostream>
#include <execution>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*256)
#define RAND_0_1 (static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX))


Particle::Particle(const glm::vec2 &pos, float radius, float damping)
    : position{pos, pos}, velocity(0.0, 0.0), force(0.0), damping(damping), radius(radius), density(0.0), lambda(0.0), deltaPosition(0.0)
{

}

Wall::Wall(const glm::vec2& pos, const glm::vec2& normal)
    : position(pos), normal(normal)
{

}

PBFSystem::PBFSystem()
    : timeStep(0.0083),
      solverIteration(2),
      damping(0.0), radius(0.075), kernelRadius(2.5f*radius),
      restDensity(1.0),
      relaxEps(600.0), xsphConstant(0.001),
      surfacePressure(true), kernelSubstep(0.25f), surfaceK(0.0001), surfaceExp(4),
      densityKernel(kernelRadius), gradientKernel(kernelRadius),
      nnSearch(PARTICLE_INIT)
{
    particles.reserve(PARTICLE_INIT);

    autotuneRestDensity();
}

Particle &PBFSystem::createParticle(const glm::vec2 &pos)
{
    assert(particles.size() < PARTICLE_INIT);

    particles.emplace_back(pos, radius, damping);
    return particles.back();
}

Particle &PBFSystem::particle(unsigned int index)
{
    assert(index < particles.size());

    return particles[index];
}

const Particle &PBFSystem::particle(unsigned int index) const
{
    assert(index < particles.size());

    return particles[index];
}

unsigned int PBFSystem::index(const Particle &p) const
{
    return std::distance(particles.data(), &p);
}

unsigned int PBFSystem::index(const Particle *p) const
{
    return std::distance(particles.data(), p);
}

void PBFSystem::wall(const glm::vec2 &pos, const glm::vec2 &normal)
{
    walls.emplace_back(pos, normal);
}

void PBFSystem::clear()
{
    particles.clear();
    walls.clear();
}

void PBFSystem::update()
{
    update(timeStep);
}

void PBFSystem::update(float dt)
{
    PerfMonitor::instance().start_frame();
    {
        profile_sample(total);

        const glm::vec2 gravity{0.0, -9.81};

        /* 1-4: predict particle positions */
        {
            profile_sample(predict);

            #pragma omp parallel for schedule(static)
            for(auto& p : particles)
            {
                p.velocity += dt * (gravity - p.damping * p.velocity);
                p.position[1] = p.position[0] + dt * p.velocity;
            }
        }

        /* 5-6: find particle neighbours */
        {
            profile_sample(nnsearch);

            nnSearch.clear();
            nnSearch.fillGrid(particles, kernelRadius);
            nnSearch.fillNeighbors(particles, kernelRadius);
        }

        {
            profile_sample(iterations);

            for(int i = 0; i < solverIteration; i++)
            {
                float stiffness = 1.0f - std::pow(1.0f - 0.95f, 1.0 / static_cast<float>(i+1));

                /* 9 - 11 calculated lambda */
                {
                    profile_sample(lambda);

                    #pragma omp parallel for schedule(static)
                    for(auto& p : particles)
                    {
                        compute_lambda(p);
                    }
                }

                /* 12 - 15 compute delta position, collision and response */
                {
                    profile_sample(deltapos);

                    #pragma omp parallel for schedule(static)
                    for(auto& p : particles)
                    {
                        compute_deltaPos(p, stiffness);
                    }

                    #pragma omp parallel for schedule(static)
                    for(auto& p : particles)
                    {
                        p.position[1] += p.deltaPosition;
                    }
                }
            }

        }

        {
            profile_sample(updatepredic);

            #pragma omp parallel for schedule(static)
            for(auto& p : particles)
            {
                p.velocity = (p.position[1] - p.position[0]) / dt;
                p.position[0] = p.position[1];
            }
        }

        {
            profile_sample(viscocity);

            #pragma omp parallel for schedule(static)
            for(auto& p : particles)
            {
                compute_viscocity(p);
            }
        }

    }
}

void PBFSystem::boundary_constraints(Particle& p, float stiffness)
{

    for(auto& wall : walls)
    {
        float pen = glm::dot(p.position[1] - wall.position, wall.normal);
        if(pen < p.radius)
        {
            auto closest = p.position[1] - ( pen - p.radius ) * wall.normal;
            glm::vec x_ij = p.position[1] - closest;

            float C = glm::dot(x_ij, wall.normal);
            if(C >= 0) C = 0.0f;

            p.deltaPosition += - stiffness * C * wall.normal;
        }
    }
}

float PBFSystem::particleDamping() const
{
    return damping;
}

void PBFSystem::particleDamping(float d)
{
    damping = d;
    for(auto& p : particles)
    {
        p.damping = damping;
    }
}

float PBFSystem::particleRadius() const
{
    return radius;
}

void PBFSystem::particleRadius(float r)
{
    radius = r;
    for(auto& p : particles)
    {
        p.radius = radius;
    }

    updateKernels();
}

void PBFSystem::autotuneRestDensity()
{
    auto sampledPositions = helper::randomPositions(radius, {-kernelRadius, -kernelRadius}, {kernelRadius, kernelRadius});

    unsigned int count = 0;
    float densityEstimate = 0.0f;
    for(const auto& p : sampledPositions)
    {
        float d = glm::length(p);
        densityEstimate += densityKernel(d, d*d);

        if(densityKernel(d, d*d) > 0)
            count++;
    }

    restDensity = densityEstimate;

    std::cout << "Particle in range: " << count << std::endl;
    std::cout << "Density Estimate: " << densityEstimate << std::endl;
}

void PBFSystem::updateKernels()
{
    kernelRadius = 2.5f * radius;
    densityKernel = kernel::std(kernelRadius);
    gradientKernel = kernel::spiky(kernelRadius);
}

void PBFSystem::compute_lambda(Particle &p)
{
    float sum_grad_C_i = 0.0f;
    glm::vec2 grad_C_i = {0.0, 0.0};

    p.density = densityKernel(0, 0);

    const auto& neighbours = nnSearch.neighbors(p, particles);
    std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
    {
        const auto d = glm::length(p.position[1] - particles[j].position[1]);

        // equation 2
        p.density += densityKernel(d, d*d);

        auto grad_j = gradientKernel.gradient(p.position[1], particles[j].position[1]) / restDensity;
        grad_C_i += grad_j;
        sum_grad_C_i += glm::dot(grad_j, grad_j);
    });
    sum_grad_C_i += glm::dot(grad_C_i, grad_C_i);

    // equation 1
    // tensily instability (clamp to non negative) or solved through artificial pressure s_corr -> in delta p computation
    auto C_i = (surfacePressure) ? p.density / restDensity - 1.0f : glm::max(p.density / restDensity - 1.0f, 0.0f);

    // equation 11
    p.lambda = (C_i != 0.0f) ? -(C_i / (sum_grad_C_i + relaxEps) ) : 0.0f;
}

void PBFSystem::compute_deltaPos(Particle &p, float stiffness)
{
    p.deltaPosition = {0.0f, 0.0f};

    const auto& neighbours = nnSearch.neighbors(p, particles);
    std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
    {
        const auto p_j = particles[j];

        float d = glm::length(p.position[1] - p_j.position[1]);
        float q = kernelSubstep * kernelRadius;
        float s_corr = -surfaceK * std::pow( densityKernel(d, d*d) / densityKernel(q, q*q), surfaceExp);
        s_corr = (surfacePressure) ? s_corr : 0.0f;

        p.deltaPosition += (p.lambda + p_j.lambda + s_corr) * gradientKernel.gradient(p.position[1], p_j.position[1]);
    });

    p.deltaPosition /= restDensity;

    boundary_constraints(p, stiffness);
}

void PBFSystem::compute_viscocity(Particle& p)
{
    glm::vec2 accum = {0.0, 0.0};

    const auto& neighbours = nnSearch.neighbors(p, particles);
    std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
    {
        auto d = glm::length(particles[j].position[1] - p.position[1]);
        auto v_ij = particles[j].velocity - p.velocity;
        accum += v_ij * densityKernel(d, d*d);
    });

    p.velocity += xsphConstant * accum;
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
