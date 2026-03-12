---
name: pt-mpo-construction
description: Guide building process tensor MPOs using ITensors.jl with PT-TEMPO algorithm for non-Markovian bath compression
---

# PT-MPO Construction

## Purpose

Guide Claude when building process tensor matrix product operators (PT-MPOs) using ITensors.jl. The PT-MPO compresses non-Markovian bath effects into a tensor network with bounded bond dimension χ.

## When to use

Invoke when writing PT-TEMPO construction code, debugging PT-MPO structure or convergence issues, or converting ITensor outputs to plain Julia arrays for downstream use.

## Libraries

- `ITensors.jl` — tensor network construction. Uses `Index`, `ITensor`, `MPS`/`MPO` types, NOT raw arrays. Bond dimensions controlled via `maxdim` and `cutoff` kwargs in SVD/truncation.
- `JLD2.jl` — serialization of output arrays
- `CairoMakie.jl` — convergence plots

## Key concept: per-bath construction

One PT-MPO is constructed **per bath** (one per chromophore site):
- Spin-boson: 1 bath
- Dimer: 2 independent baths
- FMO: 7 independent baths, each with its own spectral density parameters

Each bath j has its own PT-MPO with bond dimension χ_j.

## Inputs

- Spectral density function J(ω) with parameters per bath
- Temperature T (or β)
- Time step Δt
- Number of time steps N
- SVD truncation parameters: `maxdim` (max χ) and `cutoff` (singular value threshold)

## Outputs

Save per-bath to `data/pt_mpo/{model}_bath{j}.jld2` with keys:
- `t_steps` — time step array
- `mpo_slices` — `Vector{Matrix{ComplexF64}}`, each entry is the d²×d² superoperator for one time step. This is the plain-array export from ITensors, NOT ITensor objects.
- `chi` — achieved bond dimension
- `params` — Dict of construction parameters

## Algorithm

PT-TEMPO (Strathearn et al. 2018, Jørgensen & Pollock 2019):
1. Discretize the influence functional from the spectral density
2. Build the augmented density tensor iteratively over time steps
3. Compress via SVD at each step, controlling bond dimension with `maxdim`/`cutoff`
4. Export each time-step slice as a plain Julia matrix

## Visualization

- Bond dimension χ vs SVD cutoff threshold
- Convergence: overlay dynamics for different χ values

## Validation criteria

- Reproduce OQuPy/ACE results for identical parameters
- Bond dimension χ converges as cutoff is tightened (eventually plateaus)
- Applying PT-MPO slices to a system state must recover `classical-baselines` dynamics

## Common pitfalls

- **Index ordering:** ITensors.jl distinguishes physical indices (site) from link indices (bond). Confusing them silently produces wrong contractions. Always name indices explicitly (e.g., `Index(d, "phys_in")`, `Index(d, "phys_out")`).
- **Hilbert space dimension:** For a single qubit bath coupling, d=2. Don't accidentally set d to the full multi-site dimension D.
- **SVD cutoff tuning:** Too aggressive (large cutoff) → loses bath memory, dynamics deviate. Too loose (small cutoff) → χ explodes, computation slows. Start with cutoff=1e-8 and adjust.
- **Export boundary:** Always convert ITensor objects to plain `Matrix{ComplexF64}` before saving. Downstream skills must NOT depend on ITensors.jl.

## Pipeline context

Stage 2 of 5. Depends on: `classical-baselines` (for validation). Feeds into: `kraus-extraction` (MPO slices → Kraus operators).
