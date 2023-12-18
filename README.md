<p align="center">
  <img width="720" src="https://github.com/nikolausrauch/2D-fluids/assets/13553309/bc70e90a-3ee0-41e4-919b-a72b3801be02" />
</p>

---

This repository contains proof-of-concept C++ implementations of different *fluid dynamics solvers* in 2D.   
The goal is to provide minimalistic and straightforward implementations, emphasizing ease of understanding.   
I decided to keep simulation code mostly independent (in a single file), which means that routines are at times duplicated (e.g. SPH density estimation, eos, boundary particles, ...).

## Features

- Methods
  - [x] WCSPH: Weakly Compressible SPH
  - [x] PCISPH: Predictive-Corrective Incompressible SPH
  - [x] PBF: Position Based Fluid
  - [x] IISPH: Implicit Incompressible SPH
  - [ ] DFSPH: Divergence-Free SPH
  - [ ] FLIP: Fluid Particle in Cell
  - [ ] Stable-Fluid: Eulerian based Fluid
- Utility
  - [x] 2D-Viewer with OpenGL2 and ImGui/ImPlot ([standalone repo](https://github.com/nikolausrauch/2D-viewer))
  - [x] Minimal Performance Monitoring
  - [x] Hash-Grid Nearest Neighbor Search
  - [ ] Uniform-Grid Nearest Neighbor Search
  - [x] Generic Scene Description

![image](https://github.com/nikolausrauch/2D-fluids/assets/13553309/c800fd44-cb05-4995-9773-d94492f29b39)

## Minimal Code Example
Scene setup / description is independent of the used simulator.
**Scene** provides an interface to construct geometry (Box, Circle) that can either be a boundary or fluid body (currently no dynamic boundaries).
The simulator (**WCSPH, PCISPH, PBF, IISPH**) initializes its data (e.g. fluid and ghost particles) from the description via `void Simulation::create(const Scene&)`.

```C++
/* scene description */
Scene scene;

/*left, right, bottom boundary*/
desc.add(Scene::Box{ {-5.00f, -3.50f}, { 5.00f, -3.25f}, Scene::eType::BOUNDARY });
desc.add(Scene::Box{ {-5.00f, -3.25f}, {-4.75f,  8.00f}, Scene::eType::BOUNDARY });
desc.add(Scene::Box{ { 4.75f, -3.25f}, { 5.00f,  8.00f}, Scene::eType::BOUNDARY });

/* fluid */
desc.add(Scene::Box{ {-4.72f, -3.23f}, { 4.73f,  0.00f}, Scene::eType::FLUID_BODY });
desc.add(Scene::Circle{ {0.0f, 3.5f}, 0.75f, Scene::eType::FLUID_BODY});

/* simulation handler (Implicit Incompressible SPH) */
IISPH simulation;
simulation.timeStep = 1.0 / 240.0f;
simulation.stepPerFrame = 1;

/* load scene (initialize fluid and boundary particles) */
simulation.create(scene);

/* solve dynamics -> time += stepPerFrame*timeStep */
simulation.update();
```
![iisph_example](https://github.com/nikolausrauch/2D-fluids/assets/13553309/85f111e0-f733-4f8e-a210-6c6df730322f)

[](https://github.com/nikolausrauch/2D-fluids/assets/13553309/0545076c-d79e-4f6a-8fc3-3438ab3d22e2)