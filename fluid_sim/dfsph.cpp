#include "dfsph.h"

#include "utility/perfmonitor.h"
#include "utility/particlehelper.h"

#include <cassert>
#include <algorithm>

#include <glm/glm.hpp>

#define PARTICLE_INIT (1024*128)

DFSPH::Particle::Particle(const glm::vec2 &pos)
    : position(pos), velocity{{0.0f, 0.0f}, {0.0f, 0.0f}}, forceNoPressure{0.0, 0.0},
    alpha(0.0), kappa(0.0), density(0.0), densityPred(0.0)
{

}

DFSPH::GhostParticle::GhostParticle(const glm::vec2& pos)
    : position(pos)
{

}

DFSPH::DFSPH()
    : timeStep(0.0064f), stepsPerFrame(4),
    mass(0.006f), radiusParticle(0.05), radiusKernel(kernel::radiusMultiplier*radiusParticle),
    restDensity(1.0f), viscocityConstant(0.1f),
    densitySolverIter(3), divergenceSolverIter(3),
    kernelDensity(radiusKernel), kernelViscocity(radiusKernel),
    gravity(0.0f, -9.81f),
    fluidNNsearch(PARTICLE_INIT), boundaryNNsearch(PARTICLE_INIT),
    last_vel_index(0)
{
    fluidParticles.reserve(PARTICLE_INIT);
    boundaryParticles.reserve(PARTICLE_INIT);

    autotuneParameter();
}

void DFSPH::create(const Scene& desc)
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

    /* initialization step of DFSPH */
    stepFillNNsearch();
    stepDensityApha();
}

DFSPH::Particle& DFSPH::createParticle(const glm::vec2& pos)
{
    assert(fluidParticles.size() < PARTICLE_INIT);

    fluidParticles.emplace_back(pos);
    return fluidParticles.back();
}

void DFSPH::createParticles(const std::vector<glm::vec2>& pos)
{
    fluidParticles.insert(fluidParticles.end(), pos.begin(), pos.end());
}

DFSPH::GhostParticle& DFSPH::createGhostParticle(const glm::vec2& pos)
{
    assert(boundaryParticles.size() < PARTICLE_INIT);

    boundaryParticles.emplace_back(pos);
    return boundaryParticles.back();
}

unsigned int DFSPH::index(const Particle& p) const
{
    return std::distance(fluidParticles.data(), &p);
}

unsigned int DFSPH::index(const GhostParticle& p) const
{
    return std::distance(boundaryParticles.data(), &p);
}

void DFSPH::update()
{
    PerfMonitor::instance().start_frame();
    profile_sample(dfsph_total);

    for(int i = 0; i < stepsPerFrame; i++)
    {
        update(timeStep);
    }
}

void DFSPH::update(float dt)
{
    {
        profile_sample(dfsph_force);
        stepNoPressureForce();
    }

    {
        profile_sample(dfsph_pred_vel);
        stepPredictVelocity(dt);
    }

    {
        profile_sample(dfsph_density_solve);
        stepDensitySolver(dt);
    }

    {
        profile_sample(dfsph_update_pos);
        stepUpdatePosition(dt);
    }

    {
        profile_sample(dfsph_nnfill);
        fluidNNsearch.clear();
        fluidNNsearch.fillGrid(fluidParticles, radiusKernel);
    }

    {
        profile_sample(dfsph_density_alpha);
        stepDensityApha();
    }

    {
        profile_sample(dfsph_divergence_solve);
        stepDivergenceSolver(dt);
    }

    {
        profile_sample(dfsph_update_vel);
        stepUpdateVelocity();
    }
}

void DFSPH::stepFillNNsearch()
{
    fluidNNsearch.clear();
    fluidNNsearch.fillGrid(fluidParticles, radiusKernel);
}

void DFSPH::stepDensityApha()
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        /* nn search */
        fluidNNsearch.updateNeighbour(index(p_i), fluidParticles, radiusKernel);
        boundaryNNsearch.updateGhostNeighbour(index(p_i), p_i.position, boundaryParticles, radiusKernel);

        /* density */
        p_i.density = mass * kernelDensity(0, 0);

        glm::vec2 sum_all = {0.0, 0.0};
        float sum_mag = 0.0;

        /* fluid particles */
        const auto& neighbours = fluidNNsearch.neighbors(p_i, fluidParticles);
        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
        {
            const auto& p_j = fluidParticles[j];

            const auto d = glm::length(p_i.position - p_j.position);
            p_i.density += mass * kernelDensity(d, d*d);

            auto mWij = mass * kernelDensity.gradient(p_i.position, p_j.position);
            sum_all += mWij;
            sum_mag += glm::dot(mWij, mWij);
        });

        /* boundary particles */
        const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p_i));
        std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
        {
            const auto& p_j = boundaryParticles[j];

            const auto d = glm::length(p_i.position - p_j.position);
            p_i.density += mass * kernelDensity(d, d*d);

            auto mWij = mass * kernelDensity.gradient(p_i.position, p_j.position);
            sum_all += mWij;
        });

        /* position dependent term of velocity update in density and divergence solver */
        p_i.alpha = p_i.density / glm::max(glm::dot(sum_all, sum_all) + sum_mag, 1e-4f);
    }
}

void DFSPH::stepNoPressureForce()
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        /* gravity */
        p_i.forceNoPressure = mass * gravity;

        /* viscocity */
        const auto& neighbours = fluidNNsearch.neighbors(p_i, fluidParticles);
        std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
        {
            const auto& p_j = fluidParticles[j];

            /* viscocity */
            p_i.forceNoPressure += mass * mass * viscocityConstant * radiusKernel * 88.5f
                                   * glm::dot(p_i.velocity[0] - p_j.velocity[0], p_i.position - p_j.position)
                                   / (glm::dot(p_i.position - p_j.position, p_i.position - p_j.position) + 0.1f * radiusKernel * radiusKernel)
                                   * kernelViscocity.gradient(p_i.position, p_j.position);

        });
    }
}

void DFSPH::stepUpdatePosition(float dt)
{
    int idx = last_vel_index;
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.position = p_i.position + dt * p_i.velocity[idx];
    }
}

void DFSPH::stepUpdateVelocity()
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.velocity[0] = p_i.velocity[last_vel_index];
    }
}

void DFSPH::stepPredictVelocity(float dt)
{
    #pragma omp parallel for schedule(static)
    for(auto& p_i : fluidParticles)
    {
        p_i.velocity[1] = p_i.velocity[0] + dt / mass * p_i.forceNoPressure;
    }
}


void DFSPH::stepDensitySolver(float dt)
{
    for(int i = 1; i <= densitySolverIter; i++)
    {
        float inv_dt2 = 1.0 / (dt*dt);
        int idx = i % 2;

        /* compute density prediction and the corresponding stiffness value (pressure related) */
        #pragma omp parallel for schedule(static)
        for(auto& p_i : fluidParticles)
        {
            float delta_rho = 0.0;

            /* fluid particles */
            const auto& neighbours = fluidNNsearch.neighbors(p_i, fluidParticles);
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                delta_rho += mass * glm::dot((p_i.velocity[idx] - p_j.velocity[idx]), kernelDensity.gradient(p_i.position, p_j.position) );
            });

            /* boundary particles */
            const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                delta_rho += mass * glm::dot(p_i.velocity[idx], kernelDensity.gradient(p_i.position, p_j.position) );
            });


            /* stiffness value; only consider fluid compression */
            float pred_density = p_i.density + dt * delta_rho;
            p_i.kappa = inv_dt2 * std::max(pred_density - restDensity, 0.0f) * p_i.alpha;
        }

        /* adapt velocities */
        #pragma omp parallel for schedule(static)
        for(auto& p_i : fluidParticles)
        {
            glm::vec2 corr = {0.0, 0.0};

            /* fluid particles */
            const auto& neighbours = fluidNNsearch.neighbors(p_i, fluidParticles);
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                corr += mass * (p_i.kappa / p_i.density + p_j.kappa / p_j.density) * kernelDensity.gradient(p_i.position, p_j.position);
            });

            /* boundary particles */
            const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                corr += mass * (p_i.kappa / p_i.density) * kernelDensity.gradient(p_i.position, p_j.position);
            });


            p_i.velocity[ (idx + 1) % 2 ] = p_i.velocity[idx] - dt * corr;
        }
    }

    last_vel_index = (densitySolverIter + 1) % 2;
}

void DFSPH::stepDivergenceSolver(float dt)
{
    for(int i = 0; i < divergenceSolverIter; i++)
    {
        float inv_dt = 1.0 / dt;
        int idx = (last_vel_index + i) % 2;


        /* compute approximation of density derivative and the corresponding stiffness value for the update */
        #pragma omp parallel for schedule(static)
        for(auto& p_i : fluidParticles)
        {
            float delta_rho = 0.0;

            /* fluid particles */
            const auto& neighbours = fluidNNsearch.neighbors(p_i, fluidParticles);
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                delta_rho += mass * glm::dot( (p_i.velocity[idx] - p_j.velocity[idx]), kernelDensity.gradient(p_i.position, p_j.position) );
            });

            /* boundary particles */
            const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                delta_rho += mass * glm::dot(p_i.velocity[idx], kernelDensity.gradient(p_i.position, p_j.position) );
            });

            /* avoid negative update */
            delta_rho = std::max(delta_rho, 0.0f);

            /* skip computation if density is smaller than rest density */
            float density_pred = p_i.density + dt * delta_rho;
            if(density_pred < restDensity && p_i.density < restDensity)
            {
                delta_rho = 0.0f;
            }

            p_i.kappa = inv_dt * delta_rho * p_i.alpha;
        }


        /* adapt velocities */
        #pragma omp parallel for schedule(static)
        for(auto& p_i : fluidParticles)
        {
            glm::vec2 delta_v = {0.0, 0.0};

            /* fluid particles */
            const auto& neighbours = fluidNNsearch.neighbors(p_i, fluidParticles);
            std::for_each(neighbours.begin(), neighbours.end(), [&](auto& j)
            {
                const auto& p_j = fluidParticles[j];
                delta_v += mass * (p_i.kappa / p_i.density + p_j.kappa / p_j.density) * kernelDensity.gradient(p_i.position, p_j.position);
            });

            /* boundary particles */
            const auto& boundaryNeighbours = boundaryNNsearch.neighbors(index(p_i));
            std::for_each(boundaryNeighbours.begin(), boundaryNeighbours.end(), [&](auto& j)
            {
                const auto& p_j = boundaryParticles[j];
                delta_v += mass * (p_i.kappa / p_i.density) * kernelDensity.gradient(p_i.position, p_j.position);
            });

            p_i.velocity[ (idx + 1) % 2 ] = p_i.velocity[idx] - dt * delta_v;
        }
    }

    last_vel_index = (last_vel_index + divergenceSolverIter) % 2;
}

void DFSPH::updateRadius(float radius)
{
    radiusParticle = radius;
    radiusKernel = kernel::radiusMultiplier * radius;

    kernelDensity = kernel::std(radiusKernel);
    kernelViscocity = kernel::std(radiusKernel);

    autotuneParameter();
}

void DFSPH::autotuneParameter()
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

void DFSPH::clear()
{
    fluidParticles.clear();
    boundaryParticles.clear();

    fluidNNsearch.clear();
    boundaryNNsearch.clear();
}
