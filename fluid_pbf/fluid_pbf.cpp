#include <viewer.h>
#include <utility/perfmonitor.h>
#include <cstdlib>

#include <imgui/imgui.h>

#include "pbfsystem.h"
#include "scene.h"

int main(int argc, char** argv)
{
    Viewer viewer;
    viewer.mWindow.title = "Position-based Fluid";
    viewer.mWindow.width = 1280;
    viewer.mWindow.height = 720;
    viewer.mWindow.vsync = true;
    viewer.mWindow.mHDPI = false;
    viewer.mCamera.size /= 75.0f;
    viewer.mCamera.position = {-3.0, 0.75};

    /* Simulation Handler */
    PBFSystem sim;

    auto scene = eScene::DAMM_BREAK;
    setup_scene(sim, scene);

    viewer.onUpdate([&](Window& window, double dt)
    {
        sim.update();
    });

    viewer.onDraw([&](Window& window, double dt)
    {
        /* render points */
        viewer.mRender.pointRadius = 40 * sim.particleRadius();
        viewer.drawPoints(sim.particles.begin(), sim.particles.end(), [&](const auto& p, glm::vec2& coord, glm::vec4& color)
        {
            coord = p.position[0];
            color = {0.0, 0.0, 1.0, 1.0};
            return true;
        });

        /* render ground */
        for(const auto& wall : sim.walls)
        {
            viewer.drawBoundary(wall.position, wall.normal, 10.0f, 0.5f);
        }
    });


    viewer.onGui([&](Window& window, double dt)
    {
        ImGui::SetNextWindowPos(ImVec2{16, 16});
        ImGui::SetNextWindowSize(ImVec2{window.size().x / 3.0f, window.size().y - 160.0f * (viewer.mWindow.mHDPI ? 1.5f : 1.0f)});
        ImGui::Begin("Settings", nullptr, ImGuiWindowFlags_NoSavedSettings);
        {
            ImGui::TextColored(ImVec4(0.9, 0.6, 0.6, 1.0), "FPS: %4.2f", viewer.fps());
            ImGui::TextColored(ImVec4(0.9, 0.6, 0.6, 1.0), "Num. Particles: %llu", sim.particles.size());
            {
                auto& perf_monitor = PerfMonitor::instance();
                ImGui::TextColored(ImVec4(1.0, 0.4, 0.4, 1.0), "Total : %4.3f ms", perf_monitor.time("total"));
                ImGui::TextColored(ImVec4(0.1, 0.8, 0.0, 1.0), "\tpredict: %4.3f ms",  perf_monitor.time("predict"));
                ImGui::TextColored(ImVec4(1.0, 0.5, 0.0, 1.0), "\tnnsearch: %4.3f ms", perf_monitor.time("nnsearch"));
                ImGui::TextColored(ImVec4(0.5, 0.5, 1.0, 1.0), "\tsolver: %4.3f ms", perf_monitor.time("iterations"));
                ImGui::TextColored(ImVec4(0.5, 0.5, 1.0, 1.0), "\t\tlambda: %4.3f ms", perf_monitor.time("lambda"));
                ImGui::TextColored(ImVec4(0.5, 0.5, 1.0, 1.0), "\t\tdeltapos: %4.3f ms", perf_monitor.time("deltapos"));
                ImGui::TextColored(ImVec4(0.8, 0.4, 1.0, 1.0), "\tupdate: %4.3f ms", perf_monitor.time("updatepredic"));
                ImGui::TextColored(ImVec4(1.0, 0.5, 0.8, 1.0), "\tviscocity: %4.3f ms", perf_monitor.time("updatepredic"));
            }

            ImGui::Separator();

            ImGui::TextColored(ImVec4(0.6, 0.8, 0.6, 1.0), "PBF: ");
            {
                if(ImGui::SliderFloat("damping", &sim.damping, 0, 10))
                {
                    sim.particleDamping(sim.damping);
                }

                if(ImGui::SliderFloat("radius", &sim.radius, 0.01, 0.76 / ((viewer.mWindow.mHDPI) ? 2 : 1))) /* currrently limited due to opengl' maximum on point size */
                {
                    sim.particleRadius(sim.radius);
                    sim.autotuneRestDensity();
                }

                if(ImGui::SliderFloat("rho_0", &sim.restDensity, 1.0, 4000.0)) {}
                if(ImGui::SliderFloat("relax eps", &sim.relaxEps, 0.0, 1200.0)) {}
                if(ImGui::DragFloat("xsph viscocity", &sim.xsphConstant, 0.0001, 0.0, 0.5, "%.6f")) {}

                if(ImGui::Checkbox("surface tension", &sim.surfacePressure)) {}
                if(ImGui::SliderFloat("kernel factor", &sim.kernelSubstep, 0.01, 0.99)) {}
                if(ImGui::SliderFloat("surface K", &sim.surfaceK, 0.0, 0.1, "%.6f")) {}
                if(ImGui::SliderFloat("surface exp", &sim.surfaceExp, 0.0, 4.0)) {}
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

                if(ImGui::SliderInt("solver iter", &sim.solverIteration, 0, 50)) {}
                if(ImGui::SliderFloat("Timestep", &sim.timeStep, 0.001, 0.016, "%.6f")) {}
            }

            ImGui::Separator();

            if(ImGui::Button("reload"))
            {
                setup_scene(sim, scene);
            }
            ImGui::SameLine();
            if(ImGui::Button("autotune"))
            {
                sim.updateKernels();
                sim.autotuneRestDensity();
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

    /* callback for key events */
    viewer.onKey([&] (Window& window, Keyboard& keyboard, int key, int mod, bool press)
    {
        if(!press) return;

        /* reload scene */
        if(key == GLFW_KEY_R && press) { setup_scene(sim, scene); }

        /* add particles */
        if(key == GLFW_KEY_SPACE && press)
        {
            auto sampledPositions = helper::randomPositions(sim.particleRadius(), {-0.5, 0.0}, {0.5, 1.0});
            for(const auto& pos : sampledPositions)
            {
                sim.createParticle(pos);
            }
        }

        /* close application */
        if(key == GLFW_KEY_ESCAPE) { window.close(true); }

        /* fullscreen toggle */
        if(key == GLFW_KEY_F) { window.fullscreen( !window.fullscreen() ); }

        /* zoom */
        if(key == GLFW_KEY_W) { viewer.mCamera.zoom -= 0.25f; }
        if(key == GLFW_KEY_S) { viewer.mCamera.zoom += 0.25f; }
    });


    viewer.run();


    return EXIT_SUCCESS;
}
