#include "examplescenes.h"

#include <viewer.h>

#include <iostream>

void loadScene(eScene scene, Scene& desc)
{
    desc.clear();

    if(scene == eScene::DAMM_BREAK)
    {
        desc.add(Scene::Box{ {-5.00f, -3.50f}, { 5.00f, -3.25f}, Scene::eType::BOUNDARY });
        desc.add(Scene::Box{ {-5.00f, -3.25f}, {-4.75f,  8.00f}, Scene::eType::BOUNDARY });
        desc.add(Scene::Box{ { 4.75f, -3.25f}, { 5.00f,  8.00f}, Scene::eType::BOUNDARY });
        desc.add(Scene::Box{ {-4.75f, -3.20f}, { 0.00f,  2.00f}, Scene::eType::FLUID_BODY });
    }
    else if(scene == eScene::DAMM_COLL)
    {
        desc.add(Scene::Box{ {-5.00f, -3.50f}, { 5.00f, -3.25f}, Scene::eType::BOUNDARY });
        desc.add(Scene::Box{ {-5.00f, -3.25f}, {-4.75f,  8.00f}, Scene::eType::BOUNDARY });
        desc.add(Scene::Box{ { 4.75f, -3.25f}, { 5.00f,  8.00f}, Scene::eType::BOUNDARY });

        desc.add(Scene::Box{ {-4.75f, -3.20f}, {-2.00f,  2.00f}, Scene::eType::FLUID_BODY });
        desc.add(Scene::Box{ { 2.00f, -3.20f}, { 4.75f,  2.00f}, Scene::eType::FLUID_BODY });
    }
    else if(scene == eScene::WATER_DROP)
    {
        desc.add(Scene::Box{ {-5.00f, -3.50f}, { 5.00f, -3.25f}, Scene::eType::BOUNDARY });
        desc.add(Scene::Box{ {-5.00f, -3.25f}, {-4.75f,  8.00f}, Scene::eType::BOUNDARY });
        desc.add(Scene::Box{ { 4.75f, -3.25f}, { 5.00f,  8.00f}, Scene::eType::BOUNDARY });

        desc.add(Scene::Box{ {-4.72f, -3.23f}, { 4.73f,  0.00f}, Scene::eType::FLUID_BODY });

        desc.add(Scene::Circle{ {0.0f, 3.5f}, 1.5f, Scene::eType::FLUID_BODY});
    }
    else
    {
        std::cerr << "Example Scene not implemented!" << std::endl;
    }
}

void ImGuiSelect(eScene& selectedScene, Scene& desc, Simulation& simulation)
{
    int _scene = static_cast<int>(selectedScene);
    bool changed = ImGui::RadioButton("damm break", &_scene, 0); ImGui::SameLine();
    changed = changed || ImGui::RadioButton("damm coll", &_scene, 1); ImGui::SameLine();
    changed = changed || ImGui::RadioButton("water drop", &_scene, 2);
    if(changed)
    {
        selectedScene = static_cast<eScene>(_scene);
        loadScene(selectedScene, desc);
        simulation.create(desc);
    }
}
