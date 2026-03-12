---
name: benchmark-scaling
description: Guide running the PT-VQD benchmark ladder (spin-boson → dimer → FMO) and analyzing scaling behavior with pivot signal detection
---

# Benchmark Scaling

## Purpose

Guide Claude when running the benchmark ladder, analyzing scaling behavior, and detecting pivot signals that indicate fundamental problems with the approach.

## When to use

Invoke when running benchmark simulations across system sizes, analyzing scaling results, preparing comparison figures, or when benchmark results look unexpected.

## Libraries

- All pipeline libraries (ITensors.jl, Yao.jl, etc.) via the full pipeline
- `CairoMakie.jl` — scaling plots and comparison figures
- `JLD2.jl` — structured results storage
- `DataFrames.jl` — scaling tables (optional, can use plain Dicts)

## Benchmark ladder

Run in order. Each level validates before proceeding to the next.

### Level 1: Spin-boson (1 site, 1 bath)
- **Parameters:** Ohmic spectral density, ξ=0.1, ωc=7.5, β=5
- **Compare against:** `classical-baselines` TEMPO and T-TEDOPA, Wang et al. 2024 parameters
- **Purpose:** Validate the full pipeline end-to-end on simplest case

### Level 2: Donor-acceptor dimer (2 sites, 2 baths)
- **Parameters:** Independent Drude baths, λ=35 cm⁻¹, γ=106.18 cm⁻¹, T=300K
- **Compare against:** `classical-baselines` and trapped-ion experimental results (SimChargeTransfer2025)
- **Purpose:** Proof of concept for multi-site PT-VQD with sequential bath application

### Level 3: 7-site FMO (7 sites, 7 baths)
- **Parameters:** Adolphs & Renger 2006 Hamiltonian (exact site energies and couplings)
- **Compare against:** Classical PT methods (where feasible)
- **Purpose:** Main result. Demonstrate scaling beyond 4-site limit of path-integral Kraus methods.
- **Note:** Each chromophore has its own bath with potentially different reorganization energy. Use exact literature values.

## Metrics to record at each level

| Metric | Description |
|--------|-------------|
| Kraus rank | Number of non-negligible Kraus operators per bath |
| χ (bond dim) | PT-MPO bond dimension used |
| Qubit count | System + ancilla qubits in circuit |
| Circuit depth | Number of gate layers |
| Trace distance | Error vs classical baseline reference |
| Wall-clock time | Classical PT time vs PT-VQD time |

## Outputs

Save to `data/benchmarks/{model}_scaling.jld2` with keys:
- `metrics` — Dict with all metrics above
- `comparison_table` — summary across benchmark levels

Generate figures in `figures/`:
- Kraus rank vs system size (should stay flat ≤ χ)
- Qubit count vs system size
- Trace distance vs system size
- Wall-clock time comparison bar chart

## Validation criteria

- Spin-boson and dimer: trace distance < 1% vs classical baselines
- Kraus rank bounded by χ as system size grows (core paper claim)
- 7-site FMO: runs with ~11 qubits (mid-circuit measurement mode)

## Pivot signals — STOP and flag, don't just debug

These indicate potential fundamental problems with the approach:

1. **Kraus rank NOT bounded by χ:** If rank grows with system size despite PT-MPO compression, the core scaling claim is undermined. This is not a bug to fix — it's a research finding to report.

2. **Circuit depth too large for dimer:** If the Stinespring dilation circuit is too deep even for 2 sites, the quantum circuit approach may be impractical. Consider whether divide-and-conquer Stinespring (Azevedo2025) helps.

3. **Trotter error too large:** If the system-bath splitting error is unacceptable at feasible Δt, this is a fundamental limitation of the sequential application architecture.

When any pivot signal fires, **stop the benchmark run and report findings** rather than attempting to work around the issue. These are research-level decisions.

## Common pitfalls

- **FMO Hamiltonian values:** Use exact site energies and couplings from Adolphs & Renger 2006, Table 1. Do not approximate or round.
- **Per-bath heterogeneity:** Each FMO chromophore has different bath parameters. Don't use uniform baths.
- **"Runs" ≠ "accurate":** Always verify convergence in Δt and χ at each benchmark level before claiming results. Overlay convergence curves.
- **Apples-to-oranges timing:** When comparing wall-clock times, ensure classical and quantum simulations use the same parameters and hardware context.

## Pipeline context

Stage 5 of 5. Depends on: all previous skills (full pipeline). This is the final validation and analysis stage.
