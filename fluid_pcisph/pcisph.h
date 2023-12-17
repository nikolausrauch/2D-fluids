/*******************************************************************
 *  Predictive Corrective Incompressible SPH
 *
 *
 * author: Nikolaus Rauch
 * date: 16.12.2023
 */

#pragma once

#include "kernel.h"

#include <utility/hashgrid.h>

#include <vector>

#include <glm/glm.hpp>

struct Particle
{
    Particle(const glm::vec2& pos = glm::vec2());

    glm::vec2 position[2];
    glm::vec2 velocity[2];
    glm::vec2 forceNoPress;
    glm::vec2 forcePress;

    float density;
    float densityError;
    float pressure;
};

struct GhostParticle
{
    GhostParticle(const glm::vec2& pos = glm::vec2());
    glm::vec2 position;
};

/* member access for hashgrid */
inline const glm::vec2& position(const Particle& p) { return p.position[0]; }
inline const glm::vec2& position(const GhostParticle& p) { return p.position; }


struct Boundary
{
    Boundary(const glm::vec2& pos, const glm::vec2& norm, const glm::vec2& min, const glm::vec2& max);

    bool inside(const glm::vec2& pos);
    float penetration(const glm::vec2& pos);

    glm::vec2 position;
    glm::vec2 normal;

    glm::vec2 min;
    glm::vec2 max;
};

struct PCISPH
{
    PCISPH();

    PCISPH(const PCISPH&) = delete;
    PCISPH& operator = (const PCISPH&) = delete;

    Particle& createParticle(const glm::vec2& pos = glm::vec2());
    GhostParticle& createGhostParticle(const glm::vec2& pos = glm::vec2());
    Particle& particle(unsigned int index);
    const Particle& particle(unsigned int index) const;
    unsigned int index(const Particle& p) const;
    unsigned int index(const Particle* p) const;

    void boundary(const glm::vec2& pos, const glm::vec2& normal, const glm::vec2& min, const glm::vec2& max);

    /* Careful! this invalidates all references to points */
    void clear();

    void update(float dt);
    void update();

    /* set global values for new particles -> overwrites values from individual points! */
    double particleMass() const;
    void particleMass(float mass);

    float particleDamping() const;
    void particleDamping(float damping);

    float particleRadius() const;
    void particleRadius(float radius);

    void kernelRadius(float h);
    float kernelRadius();

    /* computes mass, particle radius, from kernel radius */
    void autotuneParams();
    void computeInitialDensity();

public:
    float timeStep;
    int numPerFrame;

    /* identical mass, damping, radii for all particles */
    float mass;
    float radiusParticle;
    float radiusKernel;

    float restDensity;
    float viscocityConstant;

    /* prediction-correction */
    int numIterations;
    float errorThresh;
    float stiffnessGradFact;

    glm::vec2 gravity;

    kernel::std densityKernel;
    kernel::spiky pressureKernel;
    kernel::std viscosityKernel;

    std::vector<Particle> particles;

    std::vector<GhostParticle> particlesBoundary;
    std::vector<Boundary> boundaries;

    HashGrid<Particle, position> nnSearch;
    HashGrid<GhostParticle, position> nnSearchBoundary;
};

namespace helper
{

std::vector<glm::vec2> randomPositions(float spacing, const glm::vec2& min, const glm::vec2& max, float jitterNoise = 0.0f);

}
