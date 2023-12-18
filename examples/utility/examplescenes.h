#pragma once

#include <fluid_sim/simulation.h>

enum eScene
{
    DAMM_BREAK,
    DAMM_COLL,
    WATER_DROP
};

void ImGuiSelect(eScene& selectedScene, Scene& desc, Simulation& simulation);
void loadScene(eScene scene, Scene& desc);
