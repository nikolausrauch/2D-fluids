/*******************************************************************
 *
 * author: Nikolaus Rauch
 * date: 13.02.2021
 *
 */

#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "wcsph.h"

enum class eScene : int
{
    DAMM_BREAK,
    DAMM_COLL,
};

void setup_scene(WCSPH& sim, eScene scene)
{
    sim.clear();

    /* setup boundaries */
    sim.boundary({0, -3.5}, {0.0, 1.0}, {-5.0, -3.5 - sim.radiusKernel}, {+5.0, -3.5});
    sim.boundary({-5.0f, 1.5}, {1.0, 0.0}, {-5.0 - sim.radiusKernel, -3.5 - sim.radiusKernel}, {-5.0, 6.0});
    sim.boundary({5.0f, 1.5}, {-1.0, 0.0}, {5.0, -3.5 - sim.radiusKernel}, {5.0 + sim.radiusKernel, 6.0});

    if(scene == eScene::DAMM_BREAK)
    {
        auto sampledPositions = helper::randomPositions(sim.particleRadius(), {-4.8, -3.2}, {0.0, 2.0});
        for(const auto& pos : sampledPositions)
        {
            sim.createParticle(pos);
        }
    }
    else if(scene == eScene::DAMM_COLL)
    {
        {
            auto sampledPositions = helper::randomPositions(sim.particleRadius(), {-4.8, -3.2}, {-2.0, 2.0});
            for(const auto& pos : sampledPositions)
            {
                sim.createParticle(pos);
            }
        }

        {
            auto sampledPositions = helper::randomPositions(sim.particleRadius(), {2.0, -3.2}, {4.8, 2.0});
            for(const auto& pos : sampledPositions)
            {
                sim.createParticle(pos);
            }
        }
    }
}

