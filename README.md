# Thinking in Optimal Control

An interactive essay series — from a falling point mass to a 6-DoF
SpaceX-style rocket landing — built as [Pluto.jl](https://plutojl.org)
notebooks and published automatically to GitHub Pages.

**Read it here:** https://jeremydwong.github.io/optimal-control-pages/

## The series

1. **The Anatomy of an Optimal Control Problem** — landing a 1-D spacecraft
   with [InfiniteOpt.jl](https://github.com/infiniteopt/InfiniteOpt.jl):
   decision variables, dynamics, boundary conditions, path constraints, and
   an objective (and what happens when the objective can't tell solutions
   apart).
2. **Nonconvexity and the Lossless Trick** — rocket engines can't throttle
   to zero; the feasible set becomes an annulus, and lossless
   convexification repairs it.
3. **Landing on Mars in 3-DoF** — the full Blackmore 2013 powered-descent
   problem: mass depletion via the log-mass trick, glide-slope cones, and
   an interactive landing you can fly with sliders.
4. *(planned)* Free final time — wrapping a line search around a convex solver.
5. *(planned)* When discretization lies — a transcription-diagnostics post-mortem.
6. *(planned)* Six degrees of freedom — attitude, gimbal cones, and
   successive convexification with trust regions.

## How it works

- Each page is a **Pluto notebook** (`*.jl` in this repo). On every push, a
  GitHub Action (`.github/workflows/ExportPluto.yaml`) runs all notebooks
  with real solvers (Ipopt, Clarabel) and publishes static HTML to the
  `gh-pages` branch.
- The **interactive demos** are self-contained JS widgets: Julia solves a
  parameter sweep at export time and embeds the solution arrays, so sliders
  respond instantly on the static page — no server.
- The **"Edit or run" button** on every page launches the notebook on
  [Binder](https://mybinder.org), giving readers a free, live Julia session
  with the actual code (boot takes a few minutes).

## Companion repository

The research code these essays are built on (including the 6-DoF SCvx
implementation and the parameter studies the war stories come from) lives at
[`blackmore-opts`](https://github.com/jeremydwong/blackmore-opts).

## Running locally

```julia
import Pluto; Pluto.run()   # then open any of the *.jl notebooks
```
