/*******************************************************************
 *
 * author: Nikolaus Rauch
 * date: 13.02.2021
 *
 * Position-Based-Dynamics Handler
 *  - stores particles, and constraints
 */

#pragma once

#include "kernel.h"

#include <utility/hashgrid.h>

#include <vector>

#include <glm/glm.hpp>

struct Particle
{
    Particle(const glm::vec2& pos = glm::vec2(), float radius = 0.05, float damping = 0.02);

    glm::vec2 position[2];
    glm::vec2 velocity;
    glm::vec2 force;

    float damping;
    float radius;

    float density;
    float lambda;
    glm::vec2 deltaPosition;
};

/* member access for hashgrid */
inline const glm::vec2& position(const Particle& p) { return p.position[1]; }

struct Wall
{
    Wall(const glm::vec2& pos, const glm::vec2& normal);

    glm::vec2 position;
    glm::vec2 normal;
};

struct PBFSystem
{
    PBFSystem();

    PBFSystem(const PBFSystem&) = delete;
    PBFSystem& operator = (const PBFSystem&) = delete;

    Particle& createParticle(const glm::vec2& pos = glm::vec2());
    Particle& particle(unsigned int index);
    const Particle& particle(unsigned int index) const;
    unsigned int index(const Particle& p) const;
    unsigned int index(const Particle* p) const;

    void wall(const glm::vec2& pos, const glm::vec2& normal);


    /* Careful! this invalidates all references to points */
    void clear();

    void update(float dt);
    void update();

    /* set global values for new particles -> overwrites values from individual points! */
    float particleDamping() const;
    void particleDamping(float damping);

    float particleRadius() const;
    void particleRadius(float radius);

    /* computes mass, particle radius, from kernel radius */
    void autotuneRestDensity();
    void updateKernels();

private:
    void compute_lambda(Particle& p);
    void compute_deltaPos(Particle& p, float stiffness);
    void compute_viscocity(Particle& p);

    void boundary_constraints(Particle &p, float stiffness);

public:
    float timeStep;

    int solverIteration;

    /* identical mass, damping, radii for all particles */
    float damping;
    float radius;
    float kernelRadius;

    float restDensity;

    float relaxEps;
    float xsphConstant;

    bool surfacePressure;
    float kernelSubstep;
    float surfaceK;
    float surfaceExp;

    kernel::std densityKernel;
    kernel::spiky gradientKernel;

    std::vector<Particle> particles;
    std::vector<Wall> walls;

    HashGrid<Particle, position> nnSearch;
};

namespace helper
{

std::vector<glm::vec2> randomPositions(float spacing, const glm::vec2& min, const glm::vec2& max, float jitterNoise = 0.0f);

}
