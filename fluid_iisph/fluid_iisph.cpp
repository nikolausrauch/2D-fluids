#include <viewer.h>
#include <utility/perfmonitor.h>
#include <cstdlib>

#include "iisph.h"
#include "scene.h"

enum class eVisualize
{
    DEFAULT,
    PRESSURE,
    DENSITY,
    VELOCITY
};

int main(int argc, char** argv)
{
    Viewer viewer;
    viewer.mWindow.title = "Implicit Incompressible SPH";
    viewer.mWindow.width = 1280;
    viewer.mWindow.height = 720;
    viewer.mWindow.vsync = true;  /* Note: call to onUpdate function depends on display refresh rate */
    viewer.mWindow.mHDPI = false; /* = true for 4k (highres) displays */
    viewer.mCamera.size /= 75.0f;
    viewer.mCamera.position = {-3.0, 0.75};

    /* Simulation Handler */
    IISPH sim;

    auto scene = eScene::DAMM_BREAK;
    setup_scene(sim, scene);

    auto visual = eVisualize::DEFAULT;
    bool show_boundary = false;
    bool pause = false;

    viewer.onUpdate([&](Window& window, double dt)
    {
        if(!pause) sim.update();
    });

    viewer.onDraw([&](Window& window, double dt)
    {
        /* render points */
        viewer.mRender.pointRadius = 40 * sim.particleRadius();
        viewer.drawPoints(sim.particles.begin(), sim.particles.end(), [&](const auto& p, glm::vec2& coord, glm::vec4& color)
        {
            coord = p.position[0];
            color = {0.0, 0.0, 1.0, 1.0};

            switch (visual)
            {
            case eVisualize::PRESSURE:
                color.r = glm::min(0.05 * p.pressure, 1.0);
                color.b = 1.0 - color.r;
                break;

            case eVisualize::DENSITY:
                color.g = (p.density / sim.restDensity);
                color.b = glm::max( (p.density / sim.restDensity) - 1.0f, 0.0f);
                break;

            case eVisualize::VELOCITY:
                color.r = 0.3 * glm::length(p.velocity[0]);
                color.g = 0.3 * glm::length(p.velocity[0]);
                color.b = 0.0;
                break;

            default: break;
            }

            return true;
        });

        if(show_boundary)
        {
            viewer.drawPoints(sim.particlesBoundary.begin(), sim.particlesBoundary.end(), [&](const auto& p, glm::vec2& coord, glm::vec4& color)
            {
                coord = p.position;
                color = {0.2, 0.2, 0.2, 0.1};

                return true;
            });
        }

        /* render ground */
        for(const auto& wall : sim.boundaries)
        {
            viewer.drawBoundary(wall.position, wall.normal, 10.0f, 0.5f);
        }
    });

    viewer.onGui([&](Window& window, double dt)
    {
        ImGui::SetNextWindowPos(ImVec2{16, 16});
        ImGui::SetNextWindowSize(ImVec2{window.size().x / 3.0f, window.size().y - 160.0f * (viewer.mWindow.mHDPI ? 1.5f : 1.0f)});
        ImGui::Begin("Settings", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings);
        {
            ImGui::TextColored(ImVec4(0.9, 0.6, 0.6, 1.0), "FPS: %4.2f", viewer.fps());
            ImGui::TextColored(ImVec4(0.9, 0.6, 0.6, 1.0), "Num. Particles: %llu", sim.particles.size());
            ImGui::TextColored(ImVec4(0.9, 0.6, 0.6, 1.0), "Num. Boundary Particles: %llu", sim.particlesBoundary.size());
            {
                auto& perf_monitor = PerfMonitor::instance();
                ImGui::TextColored(ImVec4(1.0, 0.4, 0.4, 1.0), "Total : %4.3f ms", perf_monitor.time("total"));
                ImGui::TextColored(ImVec4(0.1, 0.8, 0.0, 1.0), "\tnnfill: %4.3f ms",  perf_monitor.time("nnfill"));
                ImGui::TextColored(ImVec4(1.0, 0.5, 0.0, 1.0), "\tpred advection: %4.3f ms", perf_monitor.time("predadvection"));
                ImGui::TextColored(ImVec4(0.5, 0.5, 1.0, 1.0), "\tsolver iter: %4.3f ms", perf_monitor.time("predcorr"));
                ImGui::TextColored(ImVec4(0.5, 0.5, 1.0, 1.0), "\tintegrate: %4.3f ms", perf_monitor.time("integrate"));
            }

            ImGui::Separator();

            ImGui::TextColored(ImVec4(0.6, 0.8, 0.6, 1.0), "WSPH: ");
            {
                if(ImGui::SliderFloat("mass", &sim.mass, 0.001, 10.0))
                {
                    sim.particleMass(sim.mass);
                }

                if(ImGui::SliderFloat("radius", &sim.radiusParticle, 0.01, 2.0))
                {
                    sim.particleRadius(sim.radiusParticle);
                }

                if(ImGui::SliderFloat("kernel", &sim.radiusKernel, 0.01, 2.0))
                {
                    sim.kernelRadius(sim.radiusKernel);
                }

                ImGui::SliderFloat("rho_0", &sim.restDensity, 0.01, 100.0);
                ImGui::SliderFloat("viscocity", &sim.viscocityConstant, 0.0001, 1.0);
            }

            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.6, 0.8, 0.6, 1.0), "Scene/Numeric: ");
            {

                int _scene = static_cast<int>(scene);
                bool changed = ImGui::RadioButton("damm", &_scene, 0); ImGui::SameLine();
                changed = changed || ImGui::RadioButton("damm coll", &_scene, 1);
                if(changed)
                {
                    scene = static_cast<eScene>(_scene);
                    setup_scene(sim, scene);
                }

                ImGui::SliderFloat("dt", &sim.timeStep, 0.0001, 0.016, "%.6f");
                ImGui::SliderInt("iter/frame", &sim.numPerFrame, 1, 32);
                ImGui::SliderInt("min solver steps", &sim.numIterations, 0, 100);
            }

            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.6, 0.8, 0.6, 1.0), "Visuals: ");
            {
                int _visual = static_cast<int>(visual);
                bool changed = ImGui::RadioButton("default", &_visual, 0); ImGui::SameLine();
                changed = changed || ImGui::RadioButton("pressure", &_visual, 1); ImGui::SameLine();
                changed = changed || ImGui::RadioButton("density", &_visual, 2); ImGui::SameLine();
                changed = changed || ImGui::RadioButton("velocity", &_visual, 3);
                if(changed)
                {
                    visual = static_cast<eVisualize>(_visual);
                }

                ImGui::Checkbox("show boundary particles", &show_boundary);
            }

            if(ImGui::Button("reload"))
            {
                setup_scene(sim, scene);
            }
            ImGui::SameLine();
            if(ImGui::Button("autotune"))
            {
                sim.kernelRadius(4.0f*sim.particleRadius());
            }

        }
        ImGui::End();

        ImGui::SetNextWindowPos( ImVec2(16, window.size().y - 112 * (viewer.mWindow.mHDPI ? 1.5f : 1.0f)) );
        ImGui::SetNextWindowSize(ImVec2{window.size().x / 3.0f, -1});
        ImGui::Begin("##controls", nullptr, ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoTitleBar );
        {
            ImGui::TextColored({1.0, 1.0, 0.0, 1.0},   "[w/s keys]   ");    ImGui::SameLine(); ImGui::TextColored({1.0, 1.0, 1.0, 0.9}, " control zoom level");
            ImGui::TextColored({1.0, 1.0, 0.0, 1.0},   "[f key]      ");    ImGui::SameLine(); ImGui::TextColored({1.0, 1.0, 1.0, 0.9}, " toggle fullscreen");
            ImGui::TextColored({1.0, 1.0, 0.0, 1.0},   "[r key]      ");    ImGui::SameLine(); ImGui::TextColored({1.0, 1.0, 1.0, 0.9}, " reload scene");
            ImGui::TextColored({1.0, 1.0, 0.0, 1.0},   "[space key]  ");    ImGui::SameLine(); ImGui::TextColored({1.0, 1.0, 1.0, 0.9}, " spawn water");
        }
        ImGui::End();
    });

    viewer.onKey([&] (Window& window, Keyboard& keyboard, int key, int mod, bool press)
    {

        if(!press) return;

        if(key == GLFW_KEY_P) pause = !pause;

        /* close application */
        if(key == GLFW_KEY_ESCAPE) { window.close(true); }

        /* fullscreen toggle */
        if(key == GLFW_KEY_F) { window.fullscreen( !window.fullscreen() ); }

        /* zoom in */
        if(key == GLFW_KEY_W) { viewer.mCamera.zoom -= 0.25f; }

        /* zoom out */
        if(key == GLFW_KEY_S) { viewer.mCamera.zoom += 0.25f; }

        /* reload scene */
        if(key == GLFW_KEY_R && press) { setup_scene(sim, scene); }

        /* spawn particles */
        if(key == GLFW_KEY_SPACE && press)
        {
            auto sampledPositions = helper::randomPositions(sim.particleRadius()*2.0f, {-0.5, 0.0}, {0.5, 1.0});
            for(const auto& pos : sampledPositions)
            {
                sim.createParticle(pos);
            }
        }
    });

    viewer.run();

    return EXIT_SUCCESS;
}
