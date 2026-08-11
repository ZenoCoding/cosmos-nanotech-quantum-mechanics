<div align="center">

# COSMOS Nanotechnology · Quantum Mechanics

Computational physics projects from the California State Summer School for
Mathematics and Science (COSMOS), written in C with Python visualizations.

<img src="docs/assets/annealed-honeycomb.gif"
     alt="Optimized quantum state transfer across a honeycomb lattice"
     width="600">

</div>

## Overview

The repository contains 36 programs developed across a sequence of labs. The
early labs introduce C and numerical simulation through classical mechanics,
diffusion, random walks, matrix operations, and molecular dynamics. Lab 12
focuses on tight-binding models, quantum walks, and optimized quantum state
transfer.

| Topic | Examples |
| --- | --- |
| Numerical methods | finite differences, matrix multiplication |
| Classical mechanics | projectiles, oscillators, orbital motion |
| Many-body systems | diffusion, random walks, molecular dynamics |
| Quantum mechanics | tight-binding spectra, 1D and 2D quantum walks |
| Optimization | simulated annealing for state transfer |

## Repository layout

```text
labs/00–11/     introductory C and computational physics labs
labs/12/        quantum transport and tight-binding simulations
docs/assets/    generated plots and animations
scripts/        visualization tools
Makefile        build and verification commands
```

The original early-lab files were cumulative and have been organized into a
best-fit Lab 0–11 sequence. See [labs/README.md](labs/README.md) for the lab
index and [labs/12/README.md](labs/12/README.md) for Lab 12 details.

## Build

The C programs require a C compiler, `make`, BLAS, and LAPACK. On
Debian/Ubuntu:

```bash
sudo apt install build-essential libblas-dev liblapack-dev
make -j
```

Useful targets:

```bash
make early      # build Labs 0–11
make lab12      # build Lab 12
make check      # build all 36 programs and check the Python scripts
make clean      # remove compiled programs
```

Compiled programs are placed in `build/00/` through `build/12/`.

## Generate the visuals

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
make visuals
```

The generated plots and animations are written to `docs/assets/`.

## Run Lab 12 simulations

Run commands from the repository root after building Lab 12:

```bash
./build/12/timeevolve
./build/12/timeevolve2d 10
./build/12/timehex

./build/12/annealbarrier
./build/12/evolvebarrier

./build/12/annealhex
./build/12/evolveannealed
```

These programs update data in `labs/12/data/`. Run `make visuals` afterward to
refresh the plots and animations.
