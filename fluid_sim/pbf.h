#pragma once

#include "simulation.h"
#include "sphkernel.h"
#include "utility/hashgrid.h"

struct Particle
{
    Particle(const glm::vec2& pos = glm::vec2());

    glm::vec2 position[2];
    glm::vec2 velocity;
    glm::vec2 force;

    float density;
    float lambda;
    glm::vec2 deltaPosition;
};

struct GhostParticle
{
    GhostParticle(const glm::vec2& pos = glm::vec2());
    glm::vec2 position;
};

/* member access for hashgrid */
inline const glm::vec2& position(const Particle& p) { return p.position[1]; }
inline const glm::vec2& position(const GhostParticle& p) { return p.position; }


struct PBF : public Simulation
{
    PBF();

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
    void stepPrediction(float dt);
    void stepFillNNsearch();
    void stepComputeLambda();
    void stepDeltaPos();
    void stepUpdatePos();
    void stepIntegrate(float dt);
    void stepViscosity();

public:
    float timeStep;
    int stepsPerFrame;

    float radiusParticle;
    float radiusKernel;

    float restDensity;
    float relaxEps;
    float viscocityConstant;

    int numIterations;

    bool surfacePressure;
    float kernelSubstep;
    float surfaceK;
    float surfaceExp;

    kernel::std kernelDensity;
    kernel::spiky kernelPressure;

    glm::vec2 gravity;

    std::vector<Particle> fluidParticles;
    HashGrid<Particle, position> fluidNNsearch;

    std::vector<GhostParticle> boundaryParticles;
    HashGrid<GhostParticle, position> boundaryNNsearch;
};
