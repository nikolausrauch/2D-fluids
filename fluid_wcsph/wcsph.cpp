#include "wcsph.h"

#include <utility/perfmonitor.h>

#include <cassert>
#include <algorithm>
#include <iostream>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*256)
#define RAND_0_1 (static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX))


Particle::Particle(const glm::vec2 &pos)
    : position(pos), velocity(0.0, 0.0), force(0.0, 0.0)
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
    : timeStep(0.002), numPerFrame(5),
      mass(0.006), radiusParticle(0.05), radiusKernel(2.5f*radiusParticle),
      restDensity(1.0), eosScale(100), eosExponent(7.0), viscocityConstant(0.1f),
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

Particle& WCSPH::createBoundaryParticle(const glm::vec2& pos)
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

    auto samples = helper::randomPositions(radiusParticle / 2.0f, min, max);
    for(const auto& p : samples)
    {
        auto& particle = createBoundaryParticle(p);
        particle.density = restDensity;
    }

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
            nnSearch.updateNeighbour(p, particles, radiusKernel);
            nnSearchBoundary.updateNeighbour(std::distance(particles.data(), &p), p.position, particlesBoundary, radiusKernel);

            p.density = 0.0f;

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
                const auto& neighbours = nnSearchBoundary.neighbors(p, particles);
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

                    /* viscocity */
                    p_i.force += 2.0f * viscocityConstant * mass * mass / p_j.density * (p_i.velocity - p_j.velocity) *
                            glm::dot(p_i.position - p_j.position, viscosityKernel.gradient(p_i.position, p_j.position)) /
                            (glm::dot(p_i.position - p_j.position, p_i.position - p_j.position) + 0.01f * radiusKernel*radiusKernel );
                });
            }

            /* boundary particles */
            {
                const auto& neighbours = nnSearchBoundary.neighbors(p_i, particles);
                std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
                {
                    const auto& p_j = particlesBoundary[j];
                    p_i.force -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_i.pressure / (p_j.density*p_j.density))
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
    kernelRadius(2.5f*radiusParticle);
}

float WCSPH::pressure(float density)
{
    float p = eosScale * std::pow(density / restDensity - 1.0f, eosExponent);
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
    auto sampledPositions = helper::randomPositions(radiusParticle, {-2.0*radiusKernel, -2.0*radiusKernel}, {2.0*radiusKernel, 2.0*radiusKernel});

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

NearestNeighbor::NearestNeighbor(unsigned int maxParticles, unsigned int sizeGrid, unsigned int maxNeighbor)
{
    mHashgrid.resize( sizeGrid );
    std::for_each(mHashgrid.begin(), mHashgrid.end(), [&](auto& n){ n.reserve( maxNeighbor ); });

    mNeighbors.resize(maxParticles);
    std::for_each(mNeighbors.begin(), mNeighbors.end(), [&](auto& n){ n.reserve( maxNeighbor ); });
}

void NearestNeighbor::fillGrid(std::vector<Particle> &particles, float radius)
{
    for(unsigned int i = 0; i < particles.size(); i++)
    {
        const auto& p = particles[i];
        glm::ivec2 bucket = glm::floor(p.position / (1.0f*radius));

        unsigned int hash = ( (73856093 * bucket.x) ^ (19349663 * bucket.y) );
        mHashgrid[ hash % mHashgrid.size() ].emplace_back(i);
    }
}

void NearestNeighbor::fillNeighbors(std::vector<Particle> &particles, float radius)
{
#pragma omp parallel for schedule(static)
    for(unsigned int i = 0; i < particles.size(); i++)
    {
        updateNeighbour(i, particles, radius);
    }
}

void NearestNeighbor::updateNeighbour(unsigned int i, std::vector<Particle> &particles, float radius)
{
    const auto& p = particles[i];
    glm::ivec2 bucket = glm::floor(p.position / (1.0f*radius));

    mNeighbors[i].clear();
    for(int x = -1; x <= 1; x++)
    {
        for(int y = -1; y <= 1; y++)
        {
            glm::ivec2 search(bucket.x + x,bucket.y + y);
            unsigned int hash = ( (73856093 * search.x) ^ (19349663 * search.y) );
            for(unsigned int j : mHashgrid[ hash % mHashgrid.size() ])
            {
                if(i == j) continue;

                const auto& pj = particles[j];
                if(glm::length(p.position - pj.position) < radius)
                {
                    mNeighbors[i].emplace_back(j);
                }
            }
        }
    }
}

void NearestNeighbor::updateNeighbour(Particle &p, std::vector<Particle> &particles, float radius)
{
    updateNeighbour(std::distance(particles.data(), &p), particles, radius);
}

void NearestNeighbor::updateNeighbour(unsigned int i, const glm::vec2& pos, std::vector<Particle>& particles, float radius, bool index_ignore)
{
    glm::ivec2 bucket = glm::floor(pos / (1.0f*radius));

    mNeighbors[i].clear();
    for(int x = -1; x <= 1; x++)
    {
        for(int y = -1; y <= 1; y++)
        {
            glm::ivec2 search(bucket.x + x,bucket.y + y);
            unsigned int hash = ( (73856093 * search.x) ^ (19349663 * search.y) );
            for(unsigned int j : mHashgrid[ hash % mHashgrid.size() ])
            {
                if(!index_ignore && i == j) continue;

                const auto& pj = particles[j];
                if(glm::length(pos - pj.position) < radius)
                {
                    mNeighbors[i].emplace_back(j);
                }
            }
        }
    }
}

const std::vector<unsigned int> &NearestNeighbor::neighbors(unsigned int p_index)
{
    return mNeighbors[p_index];
}

const std::vector<unsigned int> &NearestNeighbor::neighbors(const Particle &p, const std::vector<Particle> &particles)
{
    return neighbors(std::distance(particles.data(), &p));
}

void NearestNeighbor::clear()
{
    std::for_each(mHashgrid.begin(), mHashgrid.end(), [&](auto& n){ n.clear(); });
}

unsigned int NearestNeighbor::size() const
{
    return mNeighbors.size();
}

void NearestNeighbor::resize(unsigned int size)
{
    mNeighbors.resize(size);
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
