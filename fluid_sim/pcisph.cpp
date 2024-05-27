#include "pcisph.h"

#include "utility/perfmonitor.h"
#include "utility/particlehelper.h"

#include <cassert>
#include <algorithm>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*128)

PCISPH::Particle::Particle(const glm::vec2 &pos)
    : position{pos, {0.0f, 0.0f}}, velocity{{0.0f, 0.0f}, {0.0f, 0.0f}},
      forceNoPressure(0.0f), forcePressure(0.0f),
      density(0.0f), densityError(0.0f), pressure(0.0f)
{

}

PCISPH::GhostParticle::GhostParticle(const glm::vec2& pos)
    : position(pos)
{

}


PCISPH::PCISPH()
    : timeStep(0.0032f), stepsPerFrame(8),
      mass(0.006f), radiusParticle(0.05), radiusKernel(kernel::radiusMultiplier*radiusParticle),
      restDensity(1.0f), viscocityConstant(0.1f),
      numIterations(3), errorThresh(0.01f), stiffnessGradFact(1.0f),
      kernelDensity(radiusKernel), kernelPressure(radiusKernel), kernelViscocity(radiusKernel),
      gravity(0.0f, -9.81f),
      fluidNNsearch(PARTICLE_INIT), boundaryNNsearch(PARTICLE_INIT)
{
    fluidParticles.reserve(PARTICLE_INIT);
    autotuneParameter();
}

void PCISPH::create(const Scene& desc)
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

PCISPH::Particle& PCISPH::createParticle(const glm::vec2 &pos)
{
    assert(fluidParticles.size() < PARTICLE_INIT);

    fluidParticles.emplace_back(pos);
    return fluidParticles.back();
}

void PCISPH::createParticles(const std::vector<glm::vec2>& pos)
{
    fluidParticles.insert(fluidParticles.end(), pos.begin(), pos.end());
}

PCISPH::GhostParticle& PCISPH::createGhostParticle(const glm::vec2& pos)
{
    assert(boundaryParticles.size() < PARTICLE_INIT);

    boundaryParticles.emplace_back(pos);
    return boundaryParticles.back();
}

unsigned int PCISPH::index(const Particle& p) const
{
    return std::distance(fluidParticles.data(), &p);
}

unsigned int PCISPH::index(const GhostParticle& p) const
{
    return std::distance(boundaryParticles.data(), &p);
}

void PCISPH::stepFillNNsearch()
{
    fluidNNsearch.clear();
    fluidNNsearch.fillGrid(fluidParticles, radiusKernel);
}

void PCISPH::stepNoPressureForce()
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        fluidNNsearch.updateNeighbour(index(p_i), fluidParticles, radiusKernel);
        boundaryNNsearch.updateGhostNeighbour(index(p_i), p_i.position[0], boundaryParticles, radiusKernel);

        /* gravity */
        p_i.forceNoPressure = mass * gravity;

        /* viscocity */
        const auto& neighbours = fluidNNsearch.neighbors(p_i, fluidParticles);
        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
        {
            const auto& p_j = fluidParticles[j];

            /* viscocity */
            p_i.forceNoPressure += mass * mass * viscocityConstant * radiusKernel * 88.5f
                    * glm::dot(p_i.velocity[0] - p_j.velocity[0], p_i.position[0] - p_j.position[0])
                    / (glm::dot(p_i.position[0] - p_j.position[0], p_i.position[0] - p_j.position[0]) + 0.1f * radiusKernel * radiusKernel)
                    * kernelViscocity.gradient(p_i.position[0], p_j.position[0]);

        });

        /* prepare for prediction-correction step */
        p_i.forcePressure = {0.0f, 0.0f};
        p_i.pressure = 0.0f;
    }
}

void PCISPH::stepPredict(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.velocity[1] = p_i.velocity[0] + dt / mass * (p_i.forceNoPressure + p_i.forcePressure);
        p_i.position[1] = p_i.position[0] + dt * p_i.velocity[1];
    }
}

void PCISPH::stepPressureCorrect(float k)
{
    #pragma omp parallel for schedule(static)
    for(auto& p : fluidParticles)
    {
        p.density = mass * kernelDensity(0, 0);

        /* fluid particles */
        {
            const auto& neighbours = fluidNNsearch.neighbors(index(p));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto d = glm::length(p.position[1] - fluidParticles[j].position[1]);
                p.density += mass * kernelDensity(d, d*d);
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = boundaryNNsearch.neighbors(index(p));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto d = glm::length(p.position[1] - boundaryParticles[j].position);
                p.density += mass * kernelDensity(d, d*d);
            });
        }

        p.densityError = p.density - restDensity;
        p.pressure += k * p.densityError;
        p.pressure = std::max(p.pressure, 0.0f);
    }
}

void PCISPH::stepAvgError(float& avg)
{
    avg = 0.0f;
    #pragma omp parallel for schedule(static) reduction(+:avg)
    for(auto& p : fluidParticles)
    {
        avg += std::max(p.densityError, 0.0f);
    }
    avg = avg / static_cast<float>(fluidParticles.size());
}

void PCISPH::stepPressureForce()
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.forcePressure = {0.0f, 0.0f};

        /* fluid particles */
        {
            const auto& neighbours = fluidNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                p_i.forcePressure -= mass * mass * (p_i.pressure / (p_i.density*p_i.density) + p_j.pressure / (p_j.density*p_j.density))
                        * kernelPressure.gradient(p_i.position[0], p_j.position[0]);
            });
        }

        /* boundary particles */
        {
            const auto& neighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                p_i.forcePressure -= mass * mass * (p_i.pressure / (p_i.density*p_i.density))
                        * kernelPressure.gradient(p_i.position[0], p_j.position);
            });
        }
    }
}

void PCISPH::stepIntegrate(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.velocity[0] = p_i.velocity[0] + dt / mass * (p_i.forceNoPressure + p_i.forcePressure);
        p_i.position[0] = p_i.position[0] + dt * p_i.velocity[0];
    }
}

void PCISPH::update()
{
    PerfMonitor::instance().start_frame();
    profile_sample(pcisph_total);

    for(int i = 0; i < stepsPerFrame; i++)
    {
        update(timeStep);
    }
}

void PCISPH::update(float dt)
{
    {
        profile_sample(pcisph_nnfill);
        stepFillNNsearch();
    }

    {
        profile_sample(pcisph_nopressforce);
        stepNoPressureForce();
    }

    /* correction scaling factor */
    float beta = 2.0f * (mass*mass * dt*dt / (restDensity*restDensity));
    float k = -1.0f / ( beta * -stiffnessGradFact );

    float avgError = 1.0f;
    for(int i = 0; i < numIterations && avgError > errorThresh; i++)
    {
        profile_sample(pcisph_predcorr);

        {
            profile_sample(pcisph_predict);
            stepPredict(dt);
        }

        {
            profile_sample(pcisph_presscorr);
            stepPressureCorrect(k);
        }

        stepAvgError(avgError);

        {
            profile_sample(pcisph_pressforce);
            stepPressureForce();
        }
    }

    {
        profile_sample(pcisph_integrate);
        stepIntegrate(dt);
    }
}

void PCISPH::updateRadius(float radius)
{
    radiusParticle = radius;
    radiusKernel = kernel::radiusMultiplier * radius;

    kernelDensity = kernel::std(radiusKernel);
    kernelPressure = kernel::spiky(radiusKernel);
    kernelViscocity = kernel::std(radiusKernel);

    autotuneParameter();
}

void PCISPH::autotuneParameter()
{
    auto sampledPositions = helper::randomPositions(2.0f*radiusParticle, {-2.0f*radiusKernel, -2.0f*radiusKernel}, {2.0f*radiusKernel, 2.0f*radiusKernel});

    float densityEstimate = 0.0f;
    glm::vec2 p_grad = {0.0f, 0.0f};
    glm::vec2 d_grad = {0.0f, 0.0f};
    float grad_dot_sum = 0.0f;

    for(const auto& p : sampledPositions)
    {
        float d = glm::length(p);
        densityEstimate += kernelDensity(d, d*d);

        /* exclude own contribution */
        if(glm::length2(p) < radiusParticle*radiusParticle*0.98) continue;

        /* I am unsure about this (can I actually use mixed kernel for density and pressure?) */
        auto d_wij = kernelPressure.gradient({0.0f, 0.0f}, p);
        auto p_wij = kernelPressure.gradient({0.0f, 0.0f}, p);
        p_grad += p_wij;
        d_grad += d_wij;
        grad_dot_sum += glm::dot(p_wij, d_wij);
    }

    mass = restDensity / densityEstimate;
    stiffnessGradFact = glm::dot(p_grad, d_grad) + grad_dot_sum;
}

void PCISPH::clear()
{
    fluidParticles.clear();
    boundaryParticles.clear();

    fluidNNsearch.clear();
    boundaryNNsearch.clear();
}
