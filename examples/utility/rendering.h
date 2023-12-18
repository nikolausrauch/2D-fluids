#pragma once

#include <viewer.h>
#include <vector>

enum eVisual
{
    DEFAULT = 0,
    VELOCITY,
    PRESSURE,
    DENSITY,
    FORCE
};

inline void ImGuiSelect(eVisual& selectedVisual)
{
    int _visual = static_cast<int>(selectedVisual);
    bool changed = ImGui::RadioButton("default", &_visual, 0); ImGui::SameLine();
    changed = changed || ImGui::RadioButton("velocity", &_visual, 1); ImGui::SameLine();
    changed = changed || ImGui::RadioButton("pressure", &_visual, 2); ImGui::SameLine();
    changed = changed || ImGui::RadioButton("density", &_visual, 3); ImGui::SameLine();
    changed = changed || ImGui::RadioButton("force", &_visual, 4);
    if(changed)
    {
        selectedVisual = static_cast<eVisual>(_visual);
    }
}

template<typename T>
void drawFluidParticles(Viewer& viewer, const std::vector<T>& particles, eVisual visual = eVisual::DEFAULT)
{
    viewer.drawPoints(particles.begin(), particles.end(), [&](const auto& p, glm::vec2& coord, glm::vec4& color)
    {
        coord = position(p);
        color = {0.0, 0.0, 1.0, 1.0};

        switch (visual)
        {
        case eVisual::PRESSURE:
            color.r = glm::min(0.05 * pressure(p), 1.0);
            color.b = 0.0;
            break;

        case eVisual::DENSITY:
            break;

        case eVisual::VELOCITY:
            color.r = 0.3 * glm::length(velocity(p));
            color.g = 0.3 * glm::length(velocity(p));
            color.b = 0.0;
            break;

        default: break;
        }

        return true;
    });
}

template<typename T>
void drawBoundaryParticles(Viewer& viewer, const std::vector<T>& particles)
{
    viewer.drawPoints(particles.begin(), particles.end(), [&](const auto& p, glm::vec2& coord, glm::vec4& color)
    {
        coord = position(p);
        color = {0.2, 0.2, 0.2, 0.1};

        return true;
    });
}
