#include <viewer.h>
#include <fluid_sim/wcsph.h>
#include <fluid_sim/utility/perfmonitor.h>
#include <fluid_sim/utility/particlehelper.h>

#include <utility/examplescenes.h>
#include <utility/rendering.h>

#include <cstdlib>

/* for render options in drawFluidParticles */
inline const glm::vec2& position(const WCSPH::Particle& p) { return p.position; }
inline const glm::vec2& position(const WCSPH::GhostParticle& p) { return p.position; }
inline float pressure(const WCSPH::Particle& p) { return p.pressure; }
inline const glm::vec2& velocity(const WCSPH::Particle& p) { return p.velocity; }

int main(int argc, char** argv)
{
    Viewer viewer;
    viewer.mWindow.title = "Weakly Compressible SPH";
    viewer.mWindow.width = 1280;
    viewer.mWindow.height = 720;
    viewer.mWindow.vsync = true;  /* Note: call to onUpdate function depends on display refresh rate */
    viewer.mWindow.mHDPI = false; /* = true for 4k (highres) displays */
    viewer.mCamera.size /= 75.0f;
    viewer.mCamera.position = {-3.0, 1.0};

    /* options */
    eScene selectedScene = eScene::DAMM_BREAK;
    eVisual selectedVisual = eVisual::DEFAULT;

    /* scene description */
    Scene sceneDesc;
    loadScene(selectedScene, sceneDesc);

    /* simulation handler */
    WCSPH simulation;
    simulation.create(sceneDesc);

    viewer.onUpdate([&](Window& window, double dt)
    {
        simulation.update();
    });

    viewer.onDraw([&](Window& window, double dt)
    {
        /* adjust point size; Note: does not match actual size of particles */
        viewer.mRender.pointRadius = 40 * simulation.radiusParticle;

        /* render particle data */
        drawFluidParticles(viewer, simulation.fluidParticles, selectedVisual);
        drawBoundaryParticles(viewer, simulation.boundaryParticles);
    });

    viewer.onGui([&](Window& window, double dt)
    {
        ImGui::SetNextWindowPos(ImVec2{16, 16});
        ImGui::SetNextWindowSize(ImVec2{window.size().x / 3.0f, window.size().y - 160.0f * (viewer.mWindow.mHDPI ? 1.5f : 1.0f)});
        ImGui::Begin("Settings", nullptr, ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoCollapse);
        {
            /* performance monitoring */
            ImGui::TextColored(ImVec4(0.9, 0.6, 0.6, 1.0), "FPS: %4.2f", viewer.fps());
            ImGui::TextColored(ImVec4(0.9, 0.6, 0.6, 1.0), "Num. Fluid Particles: %llu", simulation.fluidParticles.size());
            ImGui::TextColored(ImVec4(0.9, 0.6, 0.6, 1.0), "Num. Boundary Particles: %llu", simulation.boundaryParticles.size());
            {
                auto& perf_monitor = PerfMonitor::instance();
                ImGui::TextColored(ImVec4(1.0, 0.4, 0.4, 1.0), "Total : %4.3f ms", perf_monitor.time("wcsph_total"));
                ImGui::TextColored(ImVec4(0.1, 0.8, 0.0, 1.0), "\tnnfill: %4.3f ms",  perf_monitor.time("wcsph_nnfill"));
                ImGui::TextColored(ImVec4(1.0, 0.5, 0.0, 1.0), "\tdensity: %4.3f ms", perf_monitor.time("wcsph_density_pressure"));
                ImGui::TextColored(ImVec4(0.5, 0.5, 1.0, 1.0), "\tforce: %4.3f ms", perf_monitor.time("wcsph_force"));
                ImGui::TextColored(ImVec4(0.5, 0.5, 1.0, 1.0), "\tintegrate: %4.3f ms", perf_monitor.time("wcsph_integrate"));
            }

            /* SPH specific parameters */
            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.6, 0.8, 0.6, 1.0), "SPH: ");
            {
                if(ImGui::SliderFloat("radius", &simulation.radiusParticle, 0.01, 0.2))
                {
                    simulation.updateRadius(simulation.radiusParticle);
                    simulation.create(sceneDesc);
                }

                ImGui::SliderFloat("rho_0", &simulation.restDensity, 0.01, 4.0);
                ImGui::SliderFloat("eos_scale", &simulation.eosScale, 0.01, 1000);
                ImGui::SliderFloat("eos_exp", &simulation.eosExponent, 0.01, 7);
                ImGui::SliderFloat("viscocity", &simulation.viscocityConstant, 0.0001, 0.5);
            }

            /* numerical parameters */
            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.6, 0.8, 0.6, 1.0), "Numeric: ");
            {
                ImGui::SliderFloat("dt", &simulation.timeStep, 0.0001, 0.016, "%.6f");
                ImGui::SliderInt("iter/frame", &simulation.stepsPerFrame, 1, 32);
            }

            /* scene selection and render settings */
            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.6, 0.8, 0.6, 1.0), "Scene: ");
            ImGuiSelect(selectedScene, sceneDesc, simulation);

            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.6, 0.8, 0.6, 1.0), "Rendering: ");
            ImGuiSelect(selectedVisual);

        }
        ImGui::End();


        /* keyboard control information */
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
        if(!press) { return; }

        /* close application */
        if(key == GLFW_KEY_ESCAPE) { window.close(true); }

        /* fullscreen toggle */
        if(key == GLFW_KEY_F) { window.fullscreen( !window.fullscreen() ); }

        /* zoom in */
        if(key == GLFW_KEY_W) { viewer.mCamera.zoom -= 0.25f; }

        /* zoom out */
        if(key == GLFW_KEY_S) { viewer.mCamera.zoom += 0.25f; }

        /* zoom out */
        if(key == GLFW_KEY_R) { simulation.create(sceneDesc); }

        /* spawn water */
        if(key == GLFW_KEY_SPACE)
        {
            auto sampled_pos = helper::randomPositions(2.0f*simulation.radiusParticle, glm::vec2{0.0f, 3.0f}, 0.84f);
            simulation.createParticles(sampled_pos);
        }

    });

    viewer.run();

    return EXIT_SUCCESS;
}
