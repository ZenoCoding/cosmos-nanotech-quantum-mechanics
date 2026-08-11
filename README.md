<div align="center">

# Nanotechnology · Quantum Mechanics

**California State Summer School for Mathematics and Science (COSMOS)**

*A visual computational-physics laboratory in C and Python*

Thirty-six compiled programs spanning numerical methods, classical mechanics,
many-body motion, tight-binding spectra, quantum walks, and optimized quantum
state transfer.

</div>

## Quantum transport gallery

<table>
  <tr>
    <td width="50%" valign="top">
      <img src="docs/assets/quantum-chain.gif" alt="Quantum state scattering through a barrier on a chain">
      <p align="center"><strong>Barrier scattering on a 1D chain</strong></p>
      <p>This is a 1D tight-binding time evolution with a central on-site barrier. The wave function spreads from the left edge, and when it reaches the barrier some probability is reflected while some appears beyond it—a lattice analogue of quantum tunneling. This is the same broad phenomenon that our STMs rely on, although an STM tunnels through a vacuum gap rather than a single lattice site.</p>
    </td>
    <td width="50%" valign="top">
      <img src="docs/assets/annealed-barrier.gif" alt="Annealed transfer through a one-dimensional barrier">
      <p align="center"><strong>Optimized transfer through a barrier</strong></p>
      <p>If we use some special numbers for the overlaps between neighboring sites in our 1D chain, we can achieve near-perfect quantum state transfer—even through the barrier! Even though the wave function seems to spread out initially, simulated annealing chooses couplings that allow it to recombine at the far end with 98.5% probability near <em>t</em> = 8π. Engineered state transfer is analytically understood in ideal chains and may have applications in quantum communication and computing.</p>
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top">
      <img src="docs/assets/square-lattice.gif" alt="Quantum walk on a square lattice">
      <p align="center"><strong>Square-lattice quantum walk</strong></p>
      <p>Moving from a chain to a square lattice does not change the underlying quantum mechanics, but it gives the wave function many more paths to take. Starting from one central site, the probability spreads outward, interferes, and eventually reflects from the open boundaries. We use a logarithmic color scale so faint paths remain visible, with probabilities at or below 10<sup>−6</sup> shown at the same darkest color.</p>
    </td>
    <td width="50%" valign="top">
      <img src="docs/assets/honeycomb-lattice.gif" alt="Quantum walk on a honeycomb lattice">
      <p align="center"><strong>Honeycomb quantum walk</strong></p>
      <p>This is a honeycomb lattice. While slightly more complicated than a square lattice to label, it uses essentially the same tight-binding process as our 1D time-evolution problem. We build the Hamiltonian from connected neighbors and evolve the state with <em>e</em><sup>−<em>iHt</em></sup>, while the logarithmic color scale reveals faint probability spreading from the corner.</p>
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top">
      <img src="docs/assets/annealed-honeycomb.gif" alt="Optimized transfer on a honeycomb lattice">
      <p align="center"><strong>Optimized honeycomb transfer</strong></p>
      <p>We can try to create quantum state transfer on a more complicated honeycomb lattice by optimizing the overlaps between neighboring sites through simulated annealing. This probabilistic optimizer explores many coupling patterns and gets us close to a viable solution for a given lattice while keeping every original bond. In this run, it drives about 99.975% of the probability from the marked source to the marked target.</p>
    </td>
    <td width="50%" valign="top">
      <img src="docs/assets/density-of-states.png" alt="Density of states for three tight-binding models">
      <p align="center"><strong>Tight-binding density of states</strong></p>
      <p>Time evolution tells us how a state moves, while the density of states tells us which energies a lattice can support. These distributions compare a uniform chain, an alternating chain, and a square lattice. Their different shapes show how geometry and coupling patterns reorganize the available quantum spectrum.</p>
    </td>
  </tr>
</table>

## Computational physics gallery

<table>
  <tr>
    <td width="50%" valign="top">
      <img src="docs/assets/projectile-motion.gif" alt="Projectile trajectory">
      <p align="center"><strong>Projectile motion</strong></p>
      <p>This animation of a simple projectile served as our introduction to molecular dynamics. By computationally iterating through tiny time steps and calculating the force, velocity, and position updates, we discovered that we could simulate simple motion quite easily and accurately. The update is an Euler-style method—specifically first-order symplectic Euler, where position is updated before velocity—rather than the second-order Euler method from calculus.</p>
    </td>
    <td width="50%" valign="top">
      <img src="docs/assets/mass-spring.gif" alt="Mass-spring oscillator">
      <p align="center"><strong>Mass–spring oscillator</strong></p>
      <p>Our second exploration into molecular dynamics was a mass–spring system. By using Hooke’s law, <em>F</em> = −<em>kx</em>, we get the sinusoidal result of simple harmonic motion with very little effort, and the same model can include more complex effects such as drag. The displayed run is essentially undamped, which is why its amplitude remains nearly constant.</p>
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top">
      <img src="docs/assets/kepler-comparison.gif" alt="Comparison of two Kepler integration schemes">
      <p align="center"><strong>Kepler integrator comparison</strong></p>
      <p>Next, we used Newton’s law of universal gravitation to apply the same approach to a simple Kepler problem. As it turns out, the order of integration is crucial for approximately conserving energy; otherwise, quick movements can quickly make a problem unstable. The symplectic drift–kick method keeps the orbit bounded, while forward Euler rapidly drifts away.</p>
    </td>
    <td width="50%" valign="top">
      <img src="docs/assets/three-body.gif" alt="Three-body orbital dynamics">
      <p align="center"><strong>Three-body dynamics</strong></p>
      <p>With some more work, we can make all three bodies interact gravitationally with each other. The unperturbed starting condition is a special, stable, closed figure-eight solution to the three-body problem, but we add a small perturbation for a more interesting effect. The center of mass stays fixed while the trajectories slowly depart from exact repetition.</p>
    </td>
  </tr>
  <tr>
    <td width="50%" valign="top">
      <img src="docs/assets/diffusion.gif" alt="Comparison of diffusion simulations">
      <p align="center"><strong>Diffusion schemes</strong></p>
      <p>Molecular dynamics is not the only way to evolve a system one tiny step at a time: we can apply the same idea to a field using finite differences. These two runs start from the same density profile but impose fixed-value and zero-flux boundary conditions, so they separate first at the edges and then throughout the domain. The comparison shows that boundary conditions are part of the physical model, not just a detail of the implementation.</p>
    </td>
    <td width="50%" valign="top">
      <img src="docs/assets/molecular-dynamics.gif" alt="Lennard-Jones molecular dynamics">
      <p align="center"><strong>Lennard–Jones molecular dynamics</strong></p>
      <p>If we use molecular dynamics to simulate many atoms under the Lennard–Jones force, we can see the crystalline results of solid formation. As the temperature drops and the atoms lose velocity, short-range repulsion and longer-range attraction organize them into a lattice. Because the pair force is isotropic, the lowest-energy packing in two dimensions has six nearest neighbors and forms a triangular lattice.</p>
    </td>
  </tr>
</table>

## What this repository explores

The early labs build the numerical foundation: C programming, integration,
oscillators, projectile motion, diffusion, orbital dynamics, matrix
operations, three-body trajectories, and Lennard–Jones molecular dynamics.
Lab 12 applies those tools to tight-binding Hamiltonians and quantum transport.

> [!NOTE]
> The original early-lab directory was cumulative. Its files are now separated
> into a best-fit Lab 0–11 sequence based on program complexity, dependencies,
> Git history, and related UC Davis course material. See the
> [reconstruction notes](labs/README.md) for confidence levels and evidence.

| Area | Examples |
| --- | --- |
| Numerical methods | finite differences, matrix multiplication, timed computation |
| Classical dynamics | projectiles, springs, Kepler orbits, three-body motion |
| Statistical and many-body physics | random walks, diffusion, molecular dynamics |
| Quantum mechanics | tight-binding spectra, 1D and 2D quantum walks |
| Optimization | simulated annealing of couplings for state transfer |

## Repository layout

```text
.
├── Makefile                 # one-command builds and gallery generation
├── docs/assets/             # README-ready GIFs and plots
├── labs/
│   ├── 00/ … 11/            # reconstructed, topic-based early labs
│   │   ├── README.md        # topic, contents, and confidence note
│   │   ├── src/             # C programs assigned to that lab
│   │   ├── data/            # associated simulation output, when present
│   │   └── scripts/         # associated plotting tools, when present
│   └── 12/                  # original Lab 12 label
│       ├── src/             # 10 C programs plus LAPACK/BLAS helpers
│       ├── data/            # spectra, couplings, and quantum evolution
│       ├── results/         # original static figures
│       └── scripts/         # focused Lab 12 plotting tools
└── scripts/
    └── generate_gallery.py  # regenerates every visual shown above
```
Browse the [reconstructed Lab 0–11 index](labs/README.md) or the
[Lab 12 notes](labs/12/README.md).

## Build

The C programs require GCC or Clang, `make`, BLAS, and LAPACK. On
Debian/Ubuntu:

```bash
sudo apt install build-essential libblas-dev liblapack-dev
make -j
```

Successful builds are written to matching directories from `build/00/`
through `build/12/`. The build directory is intentionally ignored because
every executable is reproducible from source.

Useful targets:

```bash
make early      # compile the 26 reconstructed Lab 0–11 programs
make lab06      # compile one inferred lab by its two-digit number
make lab12      # compile the 10 Lab 12 programs
make check      # compile all 36 programs and syntax-check Python
make clean      # remove generated executables
```

## Regenerate the visuals

Create a Python environment, install the plotting dependencies, then rebuild
the full gallery:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
make visuals
```

The gallery generator samples long trajectories into compact GIFs while
leaving the full-resolution data untouched.

## Run selected Lab 12 simulations

Commands are designed to run from the repository root:

```bash
./build/12/timeevolve
./build/12/timeevolve2d 10
./build/12/timehex

./build/12/annealbarrier
./build/12/evolvebarrier

./build/12/annealhex
./build/12/evolveannealed

make visuals
```

The simulation executables update files in `labs/12/data/`; `make visuals`
then refreshes the presentation assets from those results.

## Notes

- Matrices in Lab 12 use column-major storage to interoperate with BLAS and
  LAPACK.
- Generated datasets and gallery assets are committed so the experiments are
  inspectable without rerunning expensive simulations.
- Generated machine binaries are not committed.
