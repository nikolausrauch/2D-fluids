#include "wcsph.h"

#include <utility/perfmonitor.h>

#include <cassert>
#include <algorithm>
#include <iostream>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*512)
#define RAND_0_1 (static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX))


Particle::Particle(const glm::vec2 &pos)
    : position(pos), velocity(0.0, 0.0), force(0.0, 0.0)
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


WCSPH::WCSPH()
    : timeStep(0.0016), numPerFrame(15),
      mass(0.006), radiusParticle(0.05), radiusKernel(4.0f*radiusParticle),
      restDensity(1.0), eosScale(250), eosExponent(1.0), viscocityConstant(0.15f),
      gravity(0, -9.81),
      densityKernel(radiusKernel), pressureKernel(radiusKernel), viscosityKernel(radiusKernel),
      nnSearch(PARTICLE_INIT), nnSearchBoundary(PARTICLE_INIT)
{
    particles.reserve(PARTICLE_INIT);
    autotuneMass();
}

Particle &WCSPH::createParticle(const glm::vec2 &pos)
{
    assert(particles.size() < PARTICLE_INIT);

    particles.emplace_back(pos);
    return particles.back();
}

GhostParticle& WCSPH::createGhostParticle(const glm::vec2& pos)
{
    assert(particles.size() < PARTICLE_INIT);

    particlesBoundary.emplace_back(pos);
    return particlesBoundary.back();
}

Particle &WCSPH::particle(unsigned int index)
{
    assert(index < particles.size());

    return particles[index];
}

const Particle &WCSPH::particle(unsigned int index) const
{
    assert(index < particles.size());

    return particles[index];
}

unsigned int WCSPH::index(const Particle &p) const
{
    return std::distance(particles.data(), &p);
}

unsigned int WCSPH::index(const Particle *p) const
{
    return std::distance(particles.data(), p);
}

void WCSPH::boundary(const glm::vec2 &pos, const glm::vec2 &normal, const glm::vec2& min, const glm::vec2& max)
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

void WCSPH::clear()
{
    particles.clear();
    particlesBoundary.clear();

    nnSearch.clear();
    nnSearchBoundary.clear();
}

void WCSPH::update()
{
    PerfMonitor::instance().start_frame();
    profile_sample(total);

    for(int i = 0; i < numPerFrame; i++)
    {
        update(timeStep);
    }
}

void WCSPH::update(float dt)
{
    {
        profile_sample(nnfill);
        nnSearch.clear();
        nnSearch.fillGrid(particles, radiusKernel);
    }

    /* compute densities and pressure */
    {
        profile_sample(density);

        #pragma omp parallel for schedule(static)
        for(auto& p : particles)
        {
            nnSearch.updateNeighbour(std::distance(particles.data(), &p), particles, radiusKernel);
            nnSearchBoundary.updateGhostNeighbour(std::distance(particles.data(), &p), p.position, particlesBoundary, radiusKernel);

            p.density = mass * densityKernel(0, 0);

            /* fluid particles */
            {
                const auto& neighbours = nnSearch.neighbors(p, particles);
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto d = glm::length(p.position - particles[j].position);
                    p.density += mass * densityKernel(d, d*d);
                });
            }

            /* boundary particles */
            {
                const auto& neighbours = nnSearchBoundary.neighbors(std::distance(particles.data(), &p));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto d = glm::length(p.position - particlesBoundary[j].position);
                    p.density += mass * restDensity * densityKernel(d, d*d);
                });
            }

            p.pressure = pressure(p.density);
        };
    }

    /* compute forces */
    {
        profile_sample(force);

        #pragma omp parallel for schedule(static)
        for(auto& p_i : particles)
        {
            p_i.force = mass * gravity;

            /* fluid particles */
            {
                const auto& neighbours = nnSearch.neighbors(p_i, particles);
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particles[j];

                    /* pressure gradient */
                    p_i.force -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_j.pressure / (p_j.density*p_j.density))
                            * pressureKernel.gradient(p_i.position, p_j.position);

                    const auto d = glm::length(p_i.position - p_j.position);
                    if(pressureKernel.firstDerivative(d, d*d) > 0.0f) { std::cerr << "hoi" << std::endl; }

                    /* viscocity */
                    p_i.force += 2.0f * viscocityConstant * mass * mass / p_j.density * (p_i.velocity - p_j.velocity) *
                            glm::dot(p_i.position - p_j.position, viscosityKernel.gradient(p_i.position, p_j.position)) /
                            (glm::dot(p_i.position - p_j.position, p_i.position - p_j.position) + 0.01f * radiusKernel*radiusKernel );
                });
            }

            /* boundary particles */
            {
                const auto& neighbours = nnSearchBoundary.neighbors(std::distance(particles.data(), &p_i));
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particlesBoundary[j];
                    p_i.force -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_i.pressure / (restDensity*restDensity))
                            * pressureKernel.gradient(p_i.position, p_j.position);
                });
            }
        };
    }

    /* integrate over time */
    {
        profile_sample(integrate);
        #pragma omp parallel for schedule(static)
        for(auto& p_i : particles)
        {
            p_i.velocity += dt * p_i.force / mass;
            p_i.position += dt * p_i.velocity;
        };
    }
}

double WCSPH::particleMass() const
{
    return mass;
}

void WCSPH::particleMass(float m)
{
    mass = m;
}

float WCSPH::particleRadius() const
{
    return radiusParticle;
}

void WCSPH::particleRadius(float r)
{
    radiusParticle = r;
    kernelRadius(4.0f*radiusParticle);
}

float WCSPH::pressure(float density)
{
    // [Batchelor 1967]
    float p = restDensity * eosScale / eosExponent * (std::pow(density / restDensity, eosExponent) - 1.0f);
    return (p < 0.0f) ? 0.0f : p;
}

void WCSPH::kernelRadius(float h)
{
    radiusKernel = h;

    densityKernel = kernel::std(h);
    pressureKernel = kernel::spiky(h);

    autotuneMass();
}

float WCSPH::kernelRadius()
{
    return radiusKernel;
}

void WCSPH::autotuneMass()
{
    auto sampledPositions = helper::randomPositions(radiusParticle*2.0f, {-2.0*radiusKernel, -2.0*radiusKernel}, {2.0*radiusKernel, 2.0*radiusKernel});

    unsigned int count = 0;
    float densityEstimate = 0.0f;
    for(const auto& p : sampledPositions)
    {
        float d = glm::length(p);
        densityEstimate += densityKernel(d, d*d);

        if(densityKernel(d, d*d) > 0)
            count++;
    }

    mass = restDensity / densityEstimate;

    std::cout << "Particle in range: " << count << std::endl;
    std::cout << "Density Estimate: " << densityEstimate * mass << std::endl;
    std::cout << "Autotune Mass: " << mass << std::endl;
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
