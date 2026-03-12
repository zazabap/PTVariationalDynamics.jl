---
name: classical-baselines
description: Guide setup of reference dynamics using QuantumDynamics.jl (TEMPO) and MPSDynamics.jl (T-TEDOPA) for validating PT-VQD results
---

# Classical Baselines

## Purpose

Guide Claude when setting up reference open quantum system dynamics using QuantumDynamics.jl (TEMPO) and MPSDynamics.jl (T-TEDOPA). These baselines serve as ground truth for validating every subsequent pipeline stage.

## When to use

Invoke when writing code that uses QuantumDynamics.jl or MPSDynamics.jl, or when generating reference dynamics for validation against PT-VQD results.

## Libraries

- `QuantumDynamics.jl` — TEMPO path integral method
- `MPSDynamics.jl` — T-TEDOPA chain mapping method
- `CairoMakie.jl` — publication-quality plotting
- `JLD2.jl` — data serialization

## Inputs

- System Hamiltonian (spin-boson σ_z or multi-site excitonic)
- Spectral density parameters: type (Ohmic/Drude), coupling strength, cutoff frequency
- Temperature (or inverse temperature β)

## Outputs

Save to `data/baselines/{model}_{method}.jld2` with keys:
- `t` — time array
- `observables` — Dict of observable time series (e.g., `"sigma_z"`, `"P1"`, `"P2"`)
- `params` — Dict of all simulation parameters used

## Visualization

Produce for every run (save to `figures/`):
- Population dynamics P_i(t) vs time — PNG + PDF
- Convergence overlay: curves for different Δt values on same axes
- TEMPO vs T-TEDOPA comparison overlay when both methods run

Use consistent styling: shared color palette, labeled axes with units, legends.

## Validation criteria

- Cross-validate TEMPO vs T-TEDOPA: must agree within 1% trace distance
- Reproduce literature results for spin-boson benchmark (ξ=0.1, ωc=7.5, β=5)

## Common pitfalls

- **Δt sensitivity:** Always run at least 3 time steps (e.g., Δt, Δt/2, Δt/4) and overlay results to confirm convergence.
- **TEMPO memory truncation:** The `K_max` parameter controls how many time steps of bath memory are kept. Too small → inaccurate non-Markovian effects. Start large, reduce until results change.
- **T-TEDOPA chain length:** Increase chain length until observables converge. Insufficient chain length silently gives wrong results.
- **Unit conventions:** QuantumDynamics.jl and MPSDynamics.jl may use different energy/time units. Verify by comparing a known analytic limit (e.g., Markovian regime).

## Pipeline context

This is stage 1 of 5. Outputs feed into validation for all downstream stages. The `quantum-circuit` skill depends on these baselines for end-to-end comparison.

Note: The project timeline includes a Yao.jl warm-up alongside baselines. For Yao.jl guidance, see the `quantum-circuit` skill.
