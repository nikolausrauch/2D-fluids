/*******************************************************************
 *  Weakly Compressible SPH
 *
 *
 * author: Nikolaus Rauch
 * date: 13.02.2021
 */

#pragma once

#include <vector>

#include <glm/glm.hpp>

#include "kernel.h"

struct Particle
{
    Particle(const glm::vec2& pos = glm::vec2());

    glm::vec2 position;
    glm::vec2 velocity;
    glm::vec2 force;

    float density;
    float pressure;
};

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

class NearestNeighbor
{
public:
    NearestNeighbor(unsigned int maxParticles, unsigned int sizeGrid = 512*512, unsigned int maxNeighbor = 16);

    void fillGrid(std::vector<Particle>& particles, float radius);
    void fillNeighbors(std::vector<Particle>& particles, float radius);
    void updateNeighbour(unsigned int i, std::vector<Particle> &particles, float radius);
    void updateNeighbour(Particle& p, std::vector<Particle> &particles, float radius);
    void updateNeighbour(unsigned int i, const glm::vec2& pos, std::vector<Particle> &particles, float radius, bool index_ignore = true);

    const std::vector<unsigned int>& neighbors(unsigned int p_index);
    const std::vector<unsigned int>& neighbors(const Particle& p, const std::vector<Particle> &particles);

    void clear();

    unsigned int size() const;
    void resize(unsigned int size);

protected:
    std::vector< std::vector<unsigned int> > mNeighbors;
    std::vector< std::vector<unsigned int> > mHashgrid;
};

struct WCSPH
{
    WCSPH();

    WCSPH(const WCSPH&) = delete;
    WCSPH& operator = (const WCSPH&) = delete;

    Particle& createParticle(const glm::vec2& pos = glm::vec2());
    Particle& createBoundaryParticle(const glm::vec2& pos = glm::vec2());
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

    float pressure(float density);

    void kernelRadius(float h);
    float kernelRadius();

    /* computes mass, particle radius, from kernel radius */
    void autotuneMass();

public:
    float timeStep;
    int numPerFrame;

    /* identical mass, damping, radii for all particles */
    float mass;
    float radiusParticle;
    float radiusKernel;

    float restDensity;
    float eosScale;
    float eosExponent;
    float viscocityConstant;

    glm::vec2 gravity;

    kernel::std densityKernel;
    kernel::spiky pressureKernel;
    kernel::std viscosityKernel;

    std::vector<Particle> particles;

    std::vector<Particle> particlesBoundary;
    std::vector<Boundary> boundaries;

    NearestNeighbor nnSearch;
    NearestNeighbor nnSearchBoundary;
};

namespace helper
{

std::vector<glm::vec2> randomPositions(float spacing, const glm::vec2& min, const glm::vec2& max, float jitterNoise = 0.0f);

}
