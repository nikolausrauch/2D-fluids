#pragma once

#include "simulation.h"
#include "sphkernel.h"
#include "utility/hashgrid.h"

struct DFSPH : public Simulation
{
    struct Particle
    {
        Particle(const glm::vec2& pos = glm::vec2());

        glm::vec2 position;
        glm::vec2 velocity[2];
        glm::vec2 forceNoPressure;

        float alpha;
        float kappa;
        float density;
        float densityPred;
    };

    struct GhostParticle
    {
        GhostParticle(const glm::vec2& pos = glm::vec2());
        glm::vec2 position;
    };

    /* member access for hashgrid (and rendering) */
    inline static const glm::vec2& position(const Particle& p) { return p.position; }
    inline static const glm::vec2& position(const GhostParticle& p) { return p.position; }


public:
    DFSPH();

    void create(const Scene& desc);
    Particle& createParticle(const glm::vec2& pos);
    void createParticles(const std::vector<glm::vec2>& pos);
    GhostParticle& createGhostParticle(const glm::vec2& pos);

    void clear();
    void update();
    void update(float dt);

    /* computes mass, particle radius, from kernel radius */
    void updateRadius(float radius);
    void autotuneParameter();

    unsigned int index(const Particle& p) const;
    unsigned int index(const GhostParticle& p) const;

public:
    /* simulation substeps */
    void stepFillNNsearch();
    void stepDensityApha();
private:
    void stepNoPressureForce();
    void stepUpdatePosition(float dt);
    void stepUpdateVelocity();
    void stepPredictVelocity(float dt);
    void stepDensitySolver(float dt);
    void stepDivergenceSolver(float dt);

public:
    float timeStep;
    int stepsPerFrame;

    float mass;
    float radiusParticle;
    float radiusKernel;

    float restDensity;
    float viscocityConstant;

    int densitySolverIter;
    int divergenceSolverIter;

    kernel::std kernelDensity;
    kernel::std kernelViscocity;

    glm::vec2 gravity;

    std::vector<Particle> fluidParticles;
    HashGrid<Particle, position> fluidNNsearch;

    std::vector<GhostParticle> boundaryParticles;
    HashGrid<GhostParticle, position> boundaryNNsearch;

private:
    int last_vel_index;
};
