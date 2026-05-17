# DRX_CA

A C++ implementation of **Cellular Automata–based Dynamic Recrystallization (DRX)** simulation. The simulation models grain nucleation, growth, and dislocation-density evolution on a 2D grid under a constant strain rate. A built-in SDL2 window renders the evolving microstructure live as the simulation runs — no external viewer or visualization tool required.

> Status: **Under active development.**

---

## Overview

The simulation:

1. Initializes a 2D grid of `grain_cell` objects with random orientations and dislocation densities (Monte Carlo grain initialization).
2. Iteratively applies a strain increment, updating each cell's dislocation density.
3. Once a critical strain `EPS_CR` is reached, identifies potential nucleation sites and probabilistically nucleates new grains.
4. Propagates recrystallized grain boundaries into neighboring cells based on local misorientation and mobility.
5. Pushes each rendered frame to the integrated SDL2 viewer. The final grid is also dumped to `grid.dat` for offline analysis.

Key model parameters (cell size, strain rate, temperature, mobility constants, etc.) are defined as `#define` macros at the top of `helpers.h`.

---

## Repository Layout

| File           | Purpose                                                                 |
| -------------- | ----------------------------------------------------------------------- |
| `drx.cpp`      | Main driver. Contains the simulation loop.                              |
| `drx_grid.h`   | Declarations for `grain_cell` and `_grid_` classes.                     |
| `drx_grid.cpp` | Implementation of grid initialization, neighbour search, nucleation, and growth. |
| `helpers.h`    | Simulation constants, material parameters, and helper declarations.     |
| `helpers.cpp`  | Implementations of helper functions (`gamma_l`, `mobility`, etc.).      |
| `viewer.h`     | Public interface of the SDL2 live viewer.                               |
| `viewer.cpp`   | SDL2 implementation of the live viewer.                                 |
| `CMakeLists.txt` | Build configuration.                                                   |
| `grid.dat`     | Final grid snapshot written at the end of the run.                      |

---

## Requirements

- A C++ compiler with C++11 support or newer (`g++` ≥ 5, `clang++` ≥ 3.4)
- **CMake** ≥ 3.10
- **SDL2** development headers and library (`libsdl2-dev`)
- A POSIX-like shell (Linux or macOS recommended; on Windows use WSL or MSYS2)
- A graphical desktop session — the simulation opens a window and exits on close

### Installing dependencies

**Ubuntu / Debian:**
```bash
sudo apt update
sudo apt install build-essential cmake pkg-config libsdl2-dev
```

**Fedora / RHEL:**
```bash
sudo dnf install gcc-c++ cmake pkgconf-pkg-config SDL2-devel
```

**macOS (Homebrew):**
```bash
brew install gcc cmake pkg-config sdl2
```

---

## Build

From the repository root:

```bash
mkdir build && cd build && cmake .. && make -j
```

That's it — `cmake` will locate SDL2 (using its CMake config when available, falling back to `pkg-config`) and produce `build/drx`.

Notes:
- Defaults to a Release build. Override with `cmake -DCMAKE_BUILD_TYPE=Debug ..` for a debug build, or `-DCMAKE_BUILD_TYPE=RelWithDebInfo` for an optimized build with symbols.
- `drx.cpp` `#include`s `helpers.cpp` and `drx_grid.cpp` directly. `viewer.cpp` is compiled as its own translation unit — it is **deliberately** isolated from the project's other headers because `helpers.h` defines short macros (`b`, `alpha`, …) that collide with SDL's parameter names.
- After editing any `.cpp`, `.h`, or `CMakeLists.txt`, just rerun `make -j` from `build/`.

---

## Run

From the project root:

```bash
./build/drx
```

A window titled `DRX_CA — …` opens and immediately begins animating. The terminal prints the iteration count, max grain velocity, and average dislocation density. The simulation stops when the accumulated strain reaches `EPS_FINAL` (default `1.5`) or when you close the viewer window.

### Viewer controls

(Focus the window first.)

| Key            | Action                                                            |
| -------------- | ----------------------------------------------------------------- |
| `q` / `Esc`    | Stop the simulation and close the viewer                          |
| `s`            | Save the current frame to `drx_frame_<n>.bmp`                     |
| Close button   | Stop the simulation and close the viewer                          |

The viewer uses the `tab20` qualitative colormap (matplotlib's), with grain IDs wrapping modulo 20.

---

## Configuration

All tunable parameters live in `helpers.h`. Common knobs:

| Macro              | Meaning                              | Default     |
| ------------------ | ------------------------------------ | ----------- |
| `GRID_SIZE`        | Grid edge length (cells per side)    | `200`       |
| `STATES`           | Number of initial grain orientations | `50`        |
| `ITERATIONS_MC`    | Monte Carlo grain-init iterations    | `1000000`   |
| `CELL_SIZE`        | Physical cell size (m)               | `2.0e-6`    |
| `EPS_CR`           | Critical strain for nucleation       | `0.22`      |
| `EPS_FINAL`        | Strain at which the simulation ends  | `1.5`       |
| `EPS_RATE`         | Strain rate (1/s)                    | `0.01`      |
| `T`                | Temperature (K)                      | `725.0`     |
| `DIS_DEN_MEAN`     | Initial dislocation-density mean     | `8.0e11`    |
| `BURGERS_B`        | Burgers vector magnitude (m)         | `2.56e-10`  |
| `RENDER_EVERY`     | Refresh the viewer every N iterations | `1`        |
| `RENDER_DELAY_MS`  | Extra delay (ms) per rendered frame  | `50`        |

After editing `helpers.h`, rebuild with `make -j` from the `build/` directory.

### Tuning the animation pace

By default the viewer refreshes **every iteration** and pauses **50 ms** afterwards, so you can actually see the microstructure evolve.

- **Animation too slow?** Lower `RENDER_DELAY_MS` toward `0`, and/or raise `RENDER_EVERY` (e.g. `5`) to refresh only every Nth step.
- **Animation too fast?** Raise `RENDER_DELAY_MS` (e.g. `200`) — the simulation still runs at the same physical strain rate, you just slow the wall-clock pace.

---

## Troubleshooting

- **`fatal error: SDL2/SDL.h: No such file or directory`** — install the SDL2 development package (see *Installing dependencies*).
- **`viewer: SDL_Init failed: No available video device`** — you're on a headless machine with no display. Either run on a desktop, set up X forwarding (`ssh -X`), or use a virtual framebuffer (`xvfb-run ./drx`).
- **Build error: `to_string` / `chrono` not found** — your compiler is too old; pass `-std=c++11` (or newer) explicitly.
- **Simulation appears stuck after "Maximum Grain Velocity"** — the strain is still below `EPS_CR`; only dislocation density is updating. Lower `EPS_CR` in `helpers.h` or wait.
- **Window stays blank** — the simulation is running but `RENDER_EVERY` may be very large. Reset it to `1`.

---

## License

No license has been specified yet. Treat the code as "all rights reserved" until a `LICENSE` file is added.
