#include "wcsph.h"

#include "utility/perfmonitor.h"
#include "utility/particlehelper.h"

#include <cassert>
#include <algorithm>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*128)

WCSPH::Particle::Particle(const glm::vec2 &pos)
    : position(pos), velocity(0.0, 0.0), force(0.0, 0.0),
      density(0.0f), pressure(0.0f)
{

}

WCSPH::GhostParticle::GhostParticle(const glm::vec2& pos)
    : position(pos)
{

}


WCSPH::WCSPH()
    : timeStep(0.0016f), stepsPerFrame(15),
      mass(0.006f), radiusParticle(0.05), radiusKernel(kernel::radiusMultiplier*radiusParticle),
      restDensity(1.0f), eosScale(250.0f), eosExponent(1.0f), viscocityConstant(0.15f),
      kernelDensity(radiusKernel), kernelPressure(radiusKernel), kernelViscocity(radiusKernel),
      gravity(0.0f, -9.81f),
      fluidNNsearch(PARTICLE_INIT), boundaryNNsearch(PARTICLE_INIT)
{
    fluidParticles.reserve(PARTICLE_INIT);
    boundaryParticles.reserve(PARTICLE_INIT);

    autotuneParameter();
}

void WCSPH::create(const Scene& desc)
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

WCSPH::Particle& WCSPH::createParticle(const glm::vec2 &pos)
{
    assert(fluidParticles.size() < PARTICLE_INIT);

    fluidParticles.emplace_back(pos);
    return fluidParticles.back();
}

void WCSPH::createParticles(const std::vector<glm::vec2>& pos)
{
    fluidParticles.insert(fluidParticles.end(), pos.begin(), pos.end());
}

WCSPH::GhostParticle& WCSPH::createGhostParticle(const glm::vec2& pos)
{
    assert(boundaryParticles.size() < PARTICLE_INIT);

    boundaryParticles.emplace_back(pos);
    return boundaryParticles.back();
}

unsigned int WCSPH::index(const Particle& p) const
{
    return std::distance(fluidParticles.data(), &p);
}

unsigned int WCSPH::index(const GhostParticle& p) const
{
    return std::distance(boundaryParticles.data(), &p);
}

void WCSPH::update()
{
    PerfMonitor::instance().start_frame();
    profile_sample(wcsph_total);

    for(int i = 0; i < stepsPerFrame; i++)
    {
        update(timeStep);
    }
}

void WCSPH::update(float dt)
{
    {
        profile_sample(wcsph_nnfill);
        stepFillNNsearch();
    }

    {
        profile_sample(wcsph_density_pressure);
        stepDensityPressure();
    }

    {
        profile_sample(wcsph_force);
        stepForce();
    }

    {
        profile_sample(wcsph_integrate);
        stepIntegrate(dt);
    }
}

void WCSPH::updateRadius(float radius)
{
    radiusParticle = radius;
    radiusKernel = kernel::radiusMultiplier * radius;

    kernelDensity = kernel::std(radiusKernel);
    kernelPressure = kernel::spiky(radiusKernel);
    kernelViscocity = kernel::std(radiusKernel);

    autotuneParameter();
}

void WCSPH::autotuneParameter()
{
    auto sampledPositions = helper::randomPositions(radiusParticle*2.0f, {-2.0*radiusKernel, -2.0*radiusKernel}, {2.0*radiusKernel, 2.0*radiusKernel});

    float densityEstimate = 0.0f;
    for(const auto& p : sampledPositions)
    {
        float d = glm::length(p);
        densityEstimate += kernelDensity(d, d*d);
    }

    mass = restDensity / densityEstimate;
}

void WCSPH::stepFillNNsearch()
{
    fluidNNsearch.clear();
    fluidNNsearch.fillGrid(fluidParticles, radiusKernel);
}

void WCSPH::stepDensityPressure()
{
    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        fluidNNsearch.updateNeighbour(index(p), fluidParticles, radiusKernel);
        boundaryNNsearch.updateGhostNeighbour(index(p), p.position, boundaryParticles, radiusKernel);

        p.density = mass * kernelDensity(0, 0);

        /* fluid particles */
        const auto& fluidNeighbours = fluidNNsearch.neighbors(index(p));
        std::for_each(fluidNeighbours.begin(), fluidNeighbours.end(), [&](auto& j)
        {
            const auto d = glm::length(p.position - fluidParticles[j].position);
            p.density += mass * kernelDensity(d, d*d);
        });

        const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p));
        std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
        {
            const auto d = glm::length(p.position - boundaryParticles[j].position);
            p.density += mass * kernelDensity(d, d*d);
        });

        p.pressure = pressure(p.density);
    }
}

void WCSPH::stepForce()
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.force = mass * gravity;

        /* fluid particles */
        const auto& fluidNeighbours = fluidNNsearch.neighbors(index(p_i));
        std::for_each(fluidNeighbours.begin(), fluidNeighbours.end(), [&](auto& j)
        {
            const auto& p_j = fluidParticles[j];

            /* pressure gradient */
            p_i.force -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_j.pressure / (p_j.density*p_j.density))
                    * kernelPressure.gradient(p_i.position, p_j.position);

            /* viscocity */
            p_i.force += 2.0f * viscocityConstant * mass * mass / p_j.density * (p_i.velocity - p_j.velocity) *
                    glm::dot(p_i.position - p_j.position, kernelViscocity.gradient(p_i.position, p_j.position)) /
                    (glm::dot(p_i.position - p_j.position, p_i.position - p_j.position) + 0.01f * radiusKernel*radiusKernel );
        });

        /* boundary particles */
        const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p_i));
        std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
        {
            const auto& p_j = boundaryParticles[j];
            p_i.force -= mass * mass * (p_i.pressure / (p_i.density*p_i.density))
                    * kernelPressure.gradient(p_i.position, p_j.position);
        });
    };
}

void WCSPH::stepIntegrate(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.velocity += dt * p_i.force / mass;
        p_i.position += dt * p_i.velocity;
    };
}

float WCSPH::pressure(float density)
{
    // [Batchelor 1967]
    float p = restDensity * eosScale / eosExponent * (std::pow(density / restDensity, eosExponent) - 1.0f);
    return (p < 0.0f) ? 0.0f : p;
}

void WCSPH::clear()
{
    fluidParticles.clear();
    boundaryParticles.clear();

    fluidNNsearch.clear();
    boundaryNNsearch.clear();
}
