#pragma once

#include "simulation.h"
#include "sphkernel.h"
#include "utility/hashgrid.h"

struct PCISPH : public Simulation
{
    struct Particle
    {
        Particle(const glm::vec2& pos = glm::vec2());

        glm::vec2 position[2];
        glm::vec2 velocity[2];
        glm::vec2 forceNoPressure;
        glm::vec2 forcePressure;

        float density;
        float densityError;
        float pressure;
    };

    struct GhostParticle
    {
        GhostParticle(const glm::vec2& pos = glm::vec2());
        glm::vec2 position;
    };

    /* member access for hashgrid (and rendering) */
    inline static const glm::vec2& position(const Particle& p) { return p.position[0]; }
    inline static const glm::vec2& position(const GhostParticle& p) { return p.position; }


public:
    PCISPH();

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
    void stepNoPressureForce();
    void stepPredict(float dt);
    void stepPressureCorrect(float k);
    void stepAvgError(float& avg);
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

    /* prediction-correction */
    int numIterations;
    float errorThresh;
    float stiffnessGradFact;

    kernel::std kernelDensity;
    kernel::spiky kernelPressure;
    kernel::std kernelViscocity;

    glm::vec2 gravity;

    std::vector<Particle> fluidParticles;
    HashGrid<Particle, position> fluidNNsearch;

    std::vector<GhostParticle> boundaryParticles;
    HashGrid<GhostParticle, position> boundaryNNsearch;
};
