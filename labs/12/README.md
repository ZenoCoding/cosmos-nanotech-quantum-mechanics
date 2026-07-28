# Lab 12 — Tight-Binding Quantum Transport

This lab studies spectra and time-dependent state transfer on one-dimensional,
square, and honeycomb lattices. It also uses simulated annealing to optimize
couplings for high-fidelity transfer through a barrier and across a finite
honeycomb lattice.

## Layout

- `src/` — C simulations, matrix helpers, and LAPACK/BLAS utilities
- `data/` — eigenvalues, optimized couplings, and time-evolution results
- `scripts/` — plotting and animation tools
- `results/` — original static plots

From the repository root, run `make lab12` to compile every Lab 12 program into
`build/12/`, or `make visuals` to regenerate the gallery in `docs/assets/`.
