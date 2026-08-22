# Radiation analytic verification suite

This directory contains six small, reproducible WarpX radiation verification
cases.  Every case has its own input file, metadata, run wrapper, compact
results and initial/final/comparison plots.

| Geometry | Radiation transport | Radiation-plasma coupling |
|---|---|---|
| Planar 1-D | `1d_searchlight_beer_lambert` | `1d_matter_radiation_exchange` |
| Cartesian 2-D | `2d_searchlight_beer_lambert` | `2d_matter_radiation_exchange` |
| Radial RCYL | `rcyl_radial_beer_lambert` | `rcyl_matter_radiation_exchange` |

## Validated results

All six cases pass in a clean run with the double-precision public radiation
branch at WarpX commit `8f4c0b35b`. The compact JSON and CSV files contain the
full-precision values; the principal acceptance metrics are:

| Case | Main analytic error | Energy-ledger error |
|---|---:|---:|
| 1-D Beer-Lambert | 0 | 0 |
| 2-D Beer-Lambert | 1.1e-16 | 2.0e-16 |
| RCYL Beer-Lambert | 0 | 0 |
| 1-D matter-radiation exchange | 6.45e-5 final, 4.53e-3 history | 2.0e-16 |
| 2-D matter-radiation exchange | 6.45e-5 final, 4.53e-3 history | 1.2e-14 |
| RCYL matter-radiation exchange | 4.55e-3 final/history | 4.0e-16 |

The RCYL exchange case uses a finite plasma annulus from 2.5 to 7.5 mm. This
keeps a geometry-neutral analytic problem away from the coordinate singularity
and the particle boundary. A 64/128/256-cell study gave 1.84%, 0.916% and
0.455% final radiation errors, respectively, demonstrating first-order radial
convergence; the committed result uses 256 cells.

## Benchmark families

### Searchlight / Beer-Lambert transport

A monoenergetic photon packet traverses a homogeneous absorbing medium.  The
reference solution is

```text
E(s) = E0 exp(-alpha s),
Edeposited = E0 - E(L).
```

The 2-D case uses the standard 45-degree searchlight geometry.  Searchlight
and beam tests are widely used to verify streaming-radiation solvers; examples
include HERACLES and the Quokka radiation test suite.  Adding uniform opacity
makes the final transmission an exact Beer-Lambert benchmark instead of only
a qualitative beam-shape test.

### Radiation-matter energy exchange

A spatially uniform, stationary hybrid-electron plasma begins hot while the
radiation field begins empty.  Transport and material motion are disabled, so
only grey Planck emission/absorption acts:

```text
dEr/dt = c alpha_P (a T^4 - Er)
Cv T + Er = Cv T0.
```

The reference curve is obtained by integrating this scalar ODE to tight
tolerance.  This is the same uniform radiation-matter equilibration family
introduced by Turner & Stone and reused in modern RHD verification suites.
Running it in three geometries additionally checks planar, Cartesian-area and
cylindrical-volume accounting.

## Literature

- Turner & Stone, *A Module for Radiation Hydrodynamic Calculations with
  ZEUS-2D Using Flux-Limited Diffusion*, ApJS 135, 95 (2001),
  https://doi.org/10.1086/321779
- Wibking & Krumholz, *Quokka: a code for two-moment AMR radiation
  hydrodynamics on GPUs*, MNRAS 512, 1430 (2022),
  https://doi.org/10.1093/mnras/stac439
- Gonzalez, Audit & Huynh, *HERACLES: a three-dimensional radiation
  hydrodynamics code*, A&A 464, 429 (2007),
  https://doi.org/10.1051/0004-6361:20065486
- Su & Olson, *Benchmark results for the non-equilibrium Marshak diffusion
  problem*, JQSRT 56, 337 (1996),
  https://doi.org/10.1016/0022-4073(96)84524-9

The Su-Olson problem is documented here because it is the next standard
spatial coupling benchmark.  It is not mislabeled as implemented: the present
WarpX branch does not yet provide an imposed incoming Planck boundary or the
special `Cv proportional to T^3` closure required by that semi-analytic
solution.

## Running

The default executable paths point to the public radiation-transport WarpX
worktree.  Run everything with the visualization Python environment:

```bash
/home/tomzhu0225/venvs/warpx-viz/bin/python run_all.py
```

Run one case with:

```bash
./1d_searchlight_beer_lambert/run.sh
```

Override executables or Python when needed:

```bash
WARPX_1D=/path/to/warpx.1d WARPX_2D=/path/to/warpx.2d \
WARPX_RCYL=/path/to/warpx.rcylinder PYTHON=/path/to/python \
    ./run_all.py
```

Each case writes transient plotfiles and reduced diagnostics under
`results/raw/`.  The committed deliverables are:

- `results/summary.json` -- quantitative pass/fail metrics;
- `results/profile.csv` or `results/history.csv` -- compact numerical and
  analytic data;
- `plots/initial.png`, `plots/final.png`, `plots/comparison.png` -- immediate
  visualization.
