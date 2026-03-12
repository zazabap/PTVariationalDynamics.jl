# Classical Baselines Implementation — Design Spec

**Date:** 2026-03-12
**Scope:** Step 1 of PT-VQD pipeline — spin-boson classical baselines using TEMPO and T-TEDOPA
**Skill:** `classical-baselines`

## Overview

Implement classical reference dynamics for the spin-boson model using two independent methods: TEMPO (QuantumDynamics.jl) and T-TEDOPA (MPSDynamics.jl). Cross-validate results, run parameter sweeps for physics intuition, and produce publication-quality figures. This establishes the validation foundation for all subsequent pipeline stages.

## Package Structure

```
PTVariationalDynamics.jl/
├── Project.toml                    # Package metadata + dependencies
├── src/
│   ├── PTVariationalDynamics.jl    # Module root, exports
│   ├── spectral_densities.jl      # J(ω) functions: Ohmic, Drude
│   ├── spin_boson.jl              # Spin-boson Hamiltonian + parameters
│   ├── baselines/
│   │   ├── tempo.jl               # TEMPO wrapper (QuantumDynamics.jl)
│   │   └── ttedopa.jl             # T-TEDOPA wrapper (MPSDynamics.jl)
│   └── plotting.jl                # Shared Makie plotting utilities
├── scripts/
│   ├── run_tempo_baseline.jl      # Run TEMPO for Wang et al. params
│   ├── run_ttedopa_baseline.jl    # Run T-TEDOPA for same params
│   ├── cross_validate.jl          # Compare TEMPO vs T-TEDOPA
│   └── parameter_sweep.jl         # Sweep ξ, β, ωc
├── test/
│   ├── runtests.jl
│   ├── test_spectral_densities.jl
│   └── test_spin_boson.jl
├── data/baselines/                 # JLD2 output (gitignored)
└── figures/                        # PNG + PDF output (gitignored)
```

### Key decisions

- Spectral density functions and spin-boson parameters are standalone modules reused across all pipeline stages.
- TEMPO and T-TEDOPA wrappers are thin: translate common parameter format → package-specific API, run simulation, return standardized result struct.
- Scripts are entry points that use the package, separate from `src/`.

## Data Model

### Spin-boson Hamiltonian

The system Hamiltonian is:

```
H_s = (ε/2) σ_z + (Δ/2) σ_x
```

The system-bath coupling is via σ_z:

```
H_sb = σ_z ⊗ Σ_k g_k (a_k + a_k†)
```

where the coupling constants g_k are determined by the spectral density J(ω). The initial state is |↑⟩ (σ_z eigenstate with eigenvalue +1).

### Common types

```julia
struct SpinBosonParams
    ξ::Float64       # Kondo parameter / coupling strength (dimensionless)
    ωc::Float64      # bath cutoff frequency
    β::Float64       # inverse temperature
    Δ::Float64       # tunneling (σ_x coefficient, H_s = ε/2 σ_z + Δ/2 σ_x)
    ε::Float64       # bias (σ_z coefficient)
end

struct SimulationResult
    t::Vector{Float64}
    observables::Dict{String, Vector{Float64}}  # "sigma_z", "sigma_x", etc.
    params::SpinBosonParams
    method::String       # "TEMPO" or "T-TEDOPA"
    initial_state::String  # e.g., "up" (σ_z = +1)
    metadata::Dict{String, Any}  # method-specific: Δt, K_max, chain_length, etc.
end
```

Both wrappers return `SimulationResult`. Plotting and cross-validation functions accept `SimulationResult` without knowing which method produced it.

### Spectral density interface

Using the Leggett convention (same as Wang et al. 2024):

```julia
ohmic_spectral_density(ω; ξ, ωc) = 2π * ξ * ω * exp(-ω / ωc)
```

where ξ is the dimensionless Kondo parameter. Note: some references use α = 2ξ or omit the 2π prefactor. Each wrapper must verify it matches this convention by checking the known weak-coupling (ξ → 0) Markovian limit.

Each wrapper translates this into its package-specific bath specification.

### Serialization

JLD2 format following project data conventions:
- Path: `data/baselines/{model}_{method}.jld2`
- The full `SimulationResult` struct is saved. Keys map to struct fields: `t`, `observables`, `params`, `method`, `initial_state`, `metadata`
- For compatibility with the project-wide data convention, loading code should access `result["t"]`, `result["observables"]`, `result["params"]`

## Parameters

### Primary benchmark (Wang et al. 2024)

- Spin-boson, Ohmic spectral density
- ξ = 0.1, ωc = 7.5, β = 5.0
- Δ = 1.0 (tunneling), ε = 0.0 (unbiased)
- Observable: ⟨σ_z(t)⟩ from initial state |↑⟩
- Time range: t ∈ [0, 20] (units of 1/Δ)

Note: Verify exact parameters from Wang et al. 2024 paper. If conventions differ from ideas report values, match the paper.

### Parameter sweeps

Three sweeps, varying one parameter while holding others at Wang et al. values:

| Sweep | Values | Purpose |
|-------|--------|---------|
| Coupling ξ | 0.01, 0.1, 0.5, 1.0 | Weak to strong (Kondo) regime |
| Temperature β | 0.5, 1.0, 5.0, 20.0 | Hot to cold |
| Cutoff ωc | 1.0, 5.0, 7.5, 20.0 | Slow to fast bath |

Total: 12 parameter sets per method, 24 runs total.

### Convergence checks

Run on primary benchmark for each method. Starting values (adjust based on results):
- **Δt sensitivity:** start at Δt = 0.1 (units of 1/Δ), then Δt/2 = 0.05, Δt/4 = 0.025
- **TEMPO memory:** start at K_max = 100 memory time steps, test 50, 100, 150, 200
- **T-TEDOPA chain:** start at chain_length = 60, test 30, 60, 90, 120

## Validation Criteria

### Cross-validation gate

TEMPO and T-TEDOPA must agree within 1% trace distance on the primary benchmark:
- Primary metric: trace distance between single-qubit density matrices at each time point, `D(ρ_TEMPO(t), ρ_TTEDOPA(t)) = ½ tr|ρ₁ - ρ₂|`
- Threshold: max trace distance over all t < 0.01
- Secondary check: `max(|σ_z^TEMPO(t) - σ_z^TTEDOPA(t)|)` < 0.01 (quick sanity check)
- If they disagree, debug before proceeding to sweeps

### Convergence gates (per method)

- **Δt:** last two refinements differ by < 0.1% in max|σ_z|
- **TEMPO K_max:** stable when increasing K_max by 50%
- **T-TEDOPA chain:** stable when increasing chain length by 50%

## Deliverables

### Figures (all PNG + PDF in `figures/`)

1. Primary benchmark: ⟨σ_z(t)⟩ with TEMPO and T-TEDOPA overlaid
2. Convergence plots: Δt, K_max, chain length (one per method)
3. Parameter sweep panels: 3 figures (one per swept parameter), each with 4 curves

All figures use consistent Makie.jl styling: shared color palette, labeled axes with units, legends.

### Data

All results as JLD2 in `data/baselines/`. Reproducible — scripts regenerate everything from scratch.

### Tests

- Unit: spectral density functions against known analytic values
- Unit: SpinBosonParams construction
- Unit: SimulationResult JLD2 round-trip (save → load → compare)
- Unit: cross-validation metric on synthetic known-answer data
- Integration: TEMPO wrapper smoke test (short time range)
- Integration: T-TEDOPA wrapper smoke test (short time range)

## Scope

**In scope:** Classical baselines for spin-boson (TEMPO + T-TEDOPA), cross-validation, parameter sweeps, convergence checks, figures.

**Out of scope:** Yao.jl warm-up — this is deferred to a separate task as agreed. See the `quantum-circuit` skill for Yao.jl guidance.

## Dependencies

Use latest available versions at implementation time and pin in Project.toml:

```toml
[deps]
QuantumDynamics = "..."   # Requires TEMPO support
MPSDynamics = "..."       # Requires T-TEDOPA support
CairoMakie = "..."
JLD2 = "..."
LinearAlgebra = "..."     # stdlib

[extras]
Test = "..."              # stdlib
```

## Approach

Sequential validation with shared infrastructure:
1. Set up Julia package structure and shared types
2. Implement and validate TEMPO wrapper independently
3. Implement and validate T-TEDOPA wrapper independently
4. Cross-validate on primary benchmark
5. Run parameter sweeps
6. Generate all figures
