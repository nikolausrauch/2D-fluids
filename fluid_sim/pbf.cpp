#include "pbf.h"

#include "utility/perfmonitor.h"
#include "utility/particlehelper.h"

#include <cassert>
#include <algorithm>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*128)

Particle::Particle(const glm::vec2 &pos)
    : position{pos, {0.0f, 0.0f}}, velocity(0.0f),
      force(0.0f),
      density(0.0f), lambda(0.0f), deltaPosition(0.0f, 0.0f)
{

}

GhostParticle::GhostParticle(const glm::vec2& pos)
    : position(pos)
{

}

PBF::PBF()
    : timeStep(0.0082), stepsPerFrame(3),
      radiusParticle(0.05), radiusKernel(kernel::radiusMultiplier*radiusParticle),
      restDensity(1.0f), relaxEps(600.0f), viscocityConstant(0.001f),
      numIterations(5),
      surfacePressure(false), kernelSubstep(0.25f), surfaceK(0.0001), surfaceExp(4),
      kernelDensity(radiusKernel), kernelPressure(radiusKernel),
      gravity(0.0f, -9.81f),
      fluidNNsearch(PARTICLE_INIT), boundaryNNsearch(PARTICLE_INIT)
{
    fluidParticles.reserve(PARTICLE_INIT);
    autotuneParameter();
}

void PBF::create(const Scene& desc)
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

Particle& PBF::createParticle(const glm::vec2 &pos)
{
    assert(fluidParticles.size() < PARTICLE_INIT);

    fluidParticles.emplace_back(pos);
    return fluidParticles.back();
}

void PBF::createParticles(const std::vector<glm::vec2>& pos)
{
    fluidParticles.insert(fluidParticles.end(), pos.begin(), pos.end());
}

GhostParticle& PBF::createGhostParticle(const glm::vec2& pos)
{
    assert(boundaryParticles.size() < PARTICLE_INIT);

    boundaryParticles.emplace_back(pos);
    return boundaryParticles.back();
}

unsigned int PBF::index(const Particle& p) const
{
    return std::distance(fluidParticles.data(), &p);
}

unsigned int PBF::index(const GhostParticle& p) const
{
    return std::distance(boundaryParticles.data(), &p);
}

void PBF::stepPrediction(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        p.velocity += dt * gravity;
        p.position[1] = p.position[0] + dt * p.velocity;
    }
}

void PBF::stepFillNNsearch()
{
    fluidNNsearch.clear();
    fluidNNsearch.fillGrid(fluidParticles, radiusKernel);

    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        fluidNNsearch.updateNeighbour(index(p), fluidParticles, radiusKernel);
        boundaryNNsearch.updateGhostNeighbour(index(p), p.position[1], boundaryParticles, radiusKernel);
    }
}

void PBF::stepComputeLambda()
{
    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        float sum_grad_C_i = 0.0f;
        glm::vec2 grad_C_i = {0.0, 0.0};

        p.density = kernelDensity(0, 0);

        const auto& neighbours = fluidNNsearch.neighbors(index(p));
        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
        {
            const auto d = glm::length(p.position[1] - fluidParticles[j].position[1]);

            p.density += kernelDensity(d, d*d);

            auto grad_j = kernelPressure.gradient(p.position[1], fluidParticles[j].position[1]) / restDensity;
            grad_C_i += grad_j;
            sum_grad_C_i += glm::dot(grad_j, grad_j);
        });

        const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p));
        std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
        {
            const auto d = glm::length(p.position[1] - boundaryParticles[j].position);

            p.density += kernelDensity(d, d*d);

            auto grad_j = kernelPressure.gradient(p.position[1], boundaryParticles[j].position) / restDensity;
            grad_C_i += grad_j;
            sum_grad_C_i += glm::dot(grad_j, grad_j);
        });

        sum_grad_C_i += glm::dot(grad_C_i, grad_C_i);

        // tensily instability (clamp to non negative) or solved through artificial pressure s_corr -> in delta p computation
        auto C_i = (surfacePressure) ? p.density / restDensity - 1.0f : glm::max(p.density / restDensity - 1.0f, 0.0f);
        p.lambda = (C_i != 0.0f) ? -(C_i / (sum_grad_C_i + relaxEps) ) : 0.0f;
    }
}

void PBF::stepDeltaPos()
{
    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        p.deltaPosition = {0.0f, 0.0f};

        const auto& neighbours = fluidNNsearch.neighbors(p, fluidParticles);
        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
        {
            const auto& p_j = fluidParticles[j];

            float d = glm::length(p.position[1] - p_j.position[1]);
            float q = kernelSubstep * radiusKernel;
            float s_corr = -surfaceK * std::pow( kernelDensity(d, d*d) / kernelDensity(q, q*q), surfaceExp);
            s_corr = (surfacePressure) ? s_corr : 0.0f;

            p.deltaPosition += (p.lambda + p_j.lambda + s_corr) * kernelPressure.gradient(p.position[1], p_j.position[1]);
        });

        const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p));
        std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
        {
            const auto& p_j = boundaryParticles[j];
            p.deltaPosition += (p.lambda + p.lambda) * kernelPressure.gradient(p.position[1], p_j.position);
        });

        p.deltaPosition /= restDensity;
    }
}

void PBF::stepUpdatePos()
{
    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        p.position[1] += p.deltaPosition;
    }
}

void PBF::stepIntegrate(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        p.velocity = (p.position[1] - p.position[0]) / dt;
        p.position[0] = p.position[1];
    }
}

void PBF::stepViscosity()
{
    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        glm::vec2 accum = {0.0, 0.0};

        const auto& neighbours = fluidNNsearch.neighbors(p, fluidParticles);
        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
        {
            auto d = glm::length(fluidParticles[j].position[1] - p.position[1]);
            auto v_ij = fluidParticles[j].velocity - p.velocity;
            accum += v_ij * kernelDensity(d, d*d);
        });

        p.velocity += viscocityConstant * accum;
    }
}

void PBF::update()
{
    PerfMonitor::instance().start_frame();
    profile_sample(pbf_total);

    for(int i = 0; i < stepsPerFrame; i++)
    {
        update(timeStep);
    }
}

void PBF::update(float dt)
{
    {
        profile_sample(pbf_predict);
        stepPrediction(dt);
    }

    {
        profile_sample(pbf_fillnn);
        stepFillNNsearch();
    }

    for(int i = 0; i < numIterations; i++)
    {
        profile_sample(pbf_solver);

        {
            profile_sample(pbf_lambda);
            stepComputeLambda();
        }

        {
            profile_sample(pbf_deltapos);
            stepDeltaPos();
        }

        {
            profile_sample(pbf_updatepos);
            stepUpdatePos();
        }
    }

    {
        profile_sample(pbf_integrate);
        stepIntegrate(dt);
    }

    {
        profile_sample(pbf_viscosity);
        stepViscosity();
    }
}

void PBF::updateRadius(float radius)
{
    radiusParticle = radius;
    radiusKernel = kernel::radiusMultiplier * radius;

    kernelDensity = kernel::std(radiusKernel);
    kernelPressure = kernel::spiky(radiusKernel);

    autotuneParameter();
}

void PBF::clear()
{
    fluidParticles.clear();
    boundaryParticles.clear();

    fluidNNsearch.clear();
    boundaryNNsearch.clear();
}

void PBF::autotuneParameter()
{
    auto sampledPositions = helper::randomPositions(2.0f*radiusParticle, {-2.0f*radiusKernel, -2.0f*radiusKernel}, {2.0f*radiusKernel, 2.0f*radiusKernel});
    float densityEstimate = 0.0f;

    for(const auto& p : sampledPositions)
    {
        float d = glm::length(p);
        densityEstimate += kernelDensity(d, d*d);
    }

    restDensity = densityEstimate;
}
