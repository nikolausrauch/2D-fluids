/*******************************************************************
 *
 * author: Nikolaus Rauch
 * date: 13.02.2021
 *
 */

#pragma once

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

#include "pbfsystem.h"

enum class eScene : int
{
    DAMM_BREAK = 0,
    DAMM_COLL = 1
};

void setup_scene(PBFSystem& sim, eScene scene)
{
    sim.clear();

    /* setup boundaries */
    sim.wall({0, -3.5}, {0.0, 1.0});
    sim.wall( {-5.0f, 1.5}, {1.0, 0.0});
    sim.wall({5.0f, 1.5}, {-1.0, 0.0});

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

