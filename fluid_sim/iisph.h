#pragma once

#include "simulation.h"
#include "sphkernel.h"
#include "utility/hashgrid.h"

struct Particle
{
    Particle(const glm::vec2& pos = glm::vec2());

    glm::vec2 position;
    glm::vec2 velocity[2];
    glm::vec2 force;

    float density;
    float densityAdv;
    float pressure;

    glm::vec2 d_ii;
    glm::vec2 sum_d_ij_p_j;
    float a_ii;
    float pressure_i[2];
};

struct GhostParticle
{
    GhostParticle(const glm::vec2& pos = glm::vec2());
    glm::vec2 position;
};

/* member access for hashgrid (and rendering) */
inline const glm::vec2& position(const Particle& p) { return p.position; }
inline const glm::vec2& position(const GhostParticle& p) { return p.position; }


struct IISPH : public Simulation
{
    IISPH();

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

private:
    /* simulation substeps */
    void stepFillNNsearch();
    void stepPredictVelocity(float dt);
    void stepAdvectDensity(float dt);
    void stepComputeCoeff(float dt, int idx, int idx_2);
    void stepUpdatePressure(float dt, int idx, int idx_2);
    void stepPressureForce();
    void stepIntegrate(float dt);

public:
    float timeStep;
    int stepsPerFrame;

    float mass;
    float radiusParticle;
    float radiusKernel;

    float restDensity;
    float viscocityConstant;

    /* relaxed Jacobi solver */
    int numIterations;

    kernel::std kernelDensity;
    kernel::spiky kernelPressure;
    kernel::std kernelViscocity;

    glm::vec2 gravity;

    std::vector<Particle> fluidParticles;
    HashGrid<Particle, position> fluidNNsearch;

    std::vector<GhostParticle> boundaryParticles;
    HashGrid<GhostParticle, position> boundaryNNsearch;
};
