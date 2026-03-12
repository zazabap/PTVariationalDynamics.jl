# PT-VQD Claude Code Skills — Design Spec

**Date:** 2026-03-12
**Scope:** Five stage-specific Claude Code skills for the PTVariationalDynamics.jl computational pipeline

## Overview

Five lean, domain-aware Claude Code skills aligned to the PT-VQD pipeline stages. Each skill provides: key libraries and idioms, expected inputs/outputs, validation criteria, and common pitfalls. Skills are ~50-100 lines each and reference each other for pipeline context.

## Pipeline Flow

```
classical-baselines → pt-mpo-construction → kraus-extraction → quantum-circuit → benchmark-scaling
```

Each skill corresponds to one stage. Outputs of one stage feed as inputs to the next.

## Data Format Conventions

All skills use JLD2 as the canonical serialization format. File paths follow a consistent pattern:

| Stage | Output path pattern | Keys |
|-------|-------------------|------|
| classical-baselines | `data/baselines/{model}_{method}.jld2` | `t`, `observables`, `params` |
| pt-mpo-construction | `data/pt_mpo/{model}_bath{j}.jld2` | `t_steps`, `mpo_slices` (Vector of Matrix), `chi`, `params` |
| kraus-extraction | `data/kraus/{model}_bath{j}_step{n}.jld2` | `kraus_ops` (Vector{Matrix{ComplexF64}}), `rank`, `choi_eigenvalues` |
| quantum-circuit | `data/circuit/{model}_dynamics.jld2` | `t`, `density_matrices`, `expectation_values`, `circuit_params` |
| benchmark-scaling | `data/benchmarks/{model}_scaling.jld2` | `metrics` (Dict), `comparison_table` |

## Visualization Conventions (shared across all skills)

All visualization uses Makie.jl with consistent styling:
- Dual output: PNG (quick inspection) + PDF (paper figures), saved to `figures/`
- Shared color palette and axis label conventions across all stages
- Convergence checks: overlay curves for varying control parameters (Δt, χ, chain length)
- Each skill produces stage-appropriate plots (see individual skill sections)

## Skill 1: `classical-baselines`

**Purpose:** Guide Claude when setting up reference dynamics using QuantumDynamics.jl (TEMPO) and MPSDynamics.jl (T-TEDOPA).

**Libraries:** `QuantumDynamics.jl`, `MPSDynamics.jl`, `Makie.jl`

**Input:** System Hamiltonian, spectral density parameters (Ohmic/Drude), temperature.

**Output:** Time-series of observables (`⟨σ_z(t)⟩`, site populations `P_i(t)`) saved as JLD2 files (see Data Format Conventions).

**Visualization:** Population dynamics P_i(t) vs time, convergence checks (overlay curves for different Δt or chain lengths), and TEMPO vs T-TEDOPA comparison overlays. Follows shared visualization conventions.

**Note:** The ideas report timeline (Weeks 1-2) includes a Yao.jl warm-up alongside baselines. This skill does not cover Yao.jl — that warm-up is deferred to the `quantum-circuit` skill, which should be consulted early for familiarization.

**Validation criteria:**
- Cross-validate TEMPO vs T-TEDOPA: results must agree within 1% trace distance
- Reproduce known literature results for spin-boson (ξ=0.1, ωc=7.5, β=5)

**Common pitfalls:**
- Time step convergence — always check Δt sensitivity
- Memory truncation length in TEMPO
- Chain length convergence in T-TEDOPA
- Unit conventions between the two packages

**Trigger:** When writing code that uses QuantumDynamics.jl or MPSDynamics.jl, or generating reference dynamics for validation.

## Skill 2: `pt-mpo-construction`

**Purpose:** Guide Claude when building process tensor MPOs using ITensors.jl.

**Library:** `ITensors.jl` — uses `Index`, `ITensor`, `MPS`/`MPO` types, not raw arrays. Bond dimensions controlled via `maxdim` and `cutoff` in SVD truncation.

**Input:** Spectral density function, temperature, time step Δt, number of time steps N. One PT-MPO is constructed **per bath** (one per chromophore site). For the dimer this means 2 independent baths; for FMO, 7 independent baths each with its own spectral density parameters.

**Output:** Per-bath PT-MPOs with verified bond dimension χ. Final output format is plain Julia matrices/arrays (not ITensor objects) so downstream skills don't depend on ITensors.jl. Saved as JLD2 (see Data Format Conventions).

**Visualization:** Bond dimension χ vs truncation cutoff, convergence of dynamics with increasing χ.

**Algorithm:** PT-TEMPO — iterative construction via influence functional discretization, then compression via SVD. References: Strathearn et al. 2018, Jørgensen & Pollock 2019.

**Validation criteria:**
- Reproduce OQuPy/ACE results for the same parameters
- Bond dimension χ converges with truncation cutoff
- Applying the PT-MPO to a system state recovers the same dynamics as classical baselines

**Common pitfalls:**
- Index ordering conventions in ITensors.jl (physical vs link indices)
- Forgetting to set correct Hilbert space dimension for bath discretization
- SVD cutoff too aggressive (losing bath memory) or too loose (χ blows up)

**Trigger:** When writing PT-TEMPO construction code or debugging PT-MPO structure/convergence issues.

## Skill 3: `kraus-extraction`

**Purpose:** Guide Claude when extracting Kraus operators from PT-MPO slices via Choi matrix eigendecomposition.

**Library:** `LinearAlgebra` (standard Julia) — `eigen`, `cholesky`, SVD. No heavy dependencies needed since inputs are small matrices.

**Input:** A single time-step slice of the PT-MPO (a d² × d² matrix representing the CPTP map Φ_j for bath j).

**Output:** A set of Kraus operators `{K_j^(k)}` with rank r_j ≤ χ_j, stored as `Vector{Matrix{ComplexF64}}`.

**Algorithm:** Reshape PT-MPO slice → Choi matrix C_j, eigendecompose C_j = Σ λ_k |v_k⟩⟨v_k|, reshape eigenvectors into Kraus operators K_j^(k) = √λ_k · reshape(v_k, d, d). Discard eigenvalues below threshold.

**Validation criteria:**
- CPTP completeness: Σ_k K^(k)† K^(k) = I within numerical tolerance
- Choi matrix is positive semidefinite (all eigenvalues ≥ 0)
- Kraus rank r ≤ χ (the PT-MPO bond dimension) — core claim of the paper
- Reconstructed channel Φ(ρ) = Σ_k K^(k) ρ K^(k)† matches the original PT-MPO slice

**Common pitfalls:**
- Reshaping convention (row-major vs column-major) when converting between superoperator and Choi forms
- Numerical noise producing small negative Choi eigenvalues (clamp, don't ignore)
- Julia is column-major so reshape ordering matters. Explicit formula: form Choi via `reshape(superop, d, d, d, d)` with `permutedims((1,3,2,4))`, then `reshape` back to `d² × d²`

**Visualization:** Choi eigenvalue spectra (bar plot), Kraus rank vs χ verification plots.

**Trigger:** When writing Choi decomposition code, debugging Kraus operator extraction, or verifying CPTP properties.

## Skill 4: `quantum-circuit`

**Purpose:** Guide Claude when implementing Stinespring dilation circuits and variational dynamics in Yao.jl.

**Library:** `Yao.jl` — `chain`, `put`, `control`, `measure!`, register construction, parameter management for variational circuits.

**Input:** Kraus operators `{K_j^(k)}` from the extraction stage, system Hamiltonian parameters.

**Output:** Circuit simulation results — density matrices or expectation values at each time step.

**Two sub-tasks:**

1. **Stinespring dilation:** Embed Kraus operators into a unitary on system + ancilla qubits. For each bath j: construct unitary U_j such that Tr_ancilla[U_j (ρ ⊗ |0⟩⟨0|) U_j†] = Σ_k K_j^(k) ρ K_j^(k)†. Ancilla size: ⌈log₂(r_j)⌉ qubits. **Two modes:** (a) Mid-circuit measurement and reset for ancilla reuse across baths (~11 qubits for FMO). (b) Deterministic (no-measurement) dilation with separate ancillas per bath (~35 qubits for FMO). Implement both; default to (a) but fall back to (b) if mid-circuit measurement is unavailable or introduces excessive noise.

2. **Variational ansatz (p-VQD):** Parameterized circuit V(θ_n) for system Hamiltonian evolution. Primary method: projected variational dynamics (p-VQD, Barison et al. 2021). Cost function: fidelity F = |⟨ψ(θ_{n+1})|ψ_target⟩|² between the variational state and the target after one time step. Gradients via parameter-shift rule (compatible with Yao.jl). Convergence tolerance: F > 1 - 1e-6. Fallback: direct fidelity maximization with BFGS if p-VQD is unstable.

**Validation criteria:**
- Stinespring circuit output matches direct Kraus application (trace distance < 1e-8 on statevector sim)
- Ancilla always in separable state after measurement/reset (no entanglement leakage)
- Full time-stepping loop reproduces classical baseline dynamics for spin-boson and dimer

**Common pitfalls:**
- Yao.jl uses LSB qubit ordering (qubit 1 is least significant)
- `measure!` collapses the register in-place
- Building unitaries from Kraus operators requires careful padding when r_j is not a power of 2
- Variational parameter optimization can get stuck — use multiple initial points
- Non-physical density matrices from numerical errors: check trace = 1 and eigenvalues ≥ 0 after each time step; apply trace renormalization and eigenvalue clamping if needed

**Visualization:** Circuit fidelity vs time step, variational parameter convergence, ancilla measurement statistics.

**Dependencies:** `kraus-extraction` (for Kraus operator inputs), `classical-baselines` (for validation targets).

**Trigger:** When writing Yao.jl circuit code, implementing Stinespring dilation, or debugging variational optimization.

## Skill 5: `benchmark-scaling`

**Purpose:** Guide Claude when running the benchmark ladder and analyzing scaling behavior.

**Input:** The full PT-VQD pipeline from previous stages, system parameters for each benchmark level.

**Benchmark ladder:**
1. **Spin-boson:** Ohmic, ξ=0.1, ωc=7.5, β=5. Validate against classical baselines and Wang et al. 2024.
2. **Donor-acceptor dimer:** 2 sites, independent Drude baths, λ=35 cm⁻¹, γ=106.18 cm⁻¹, T=300K. Compare against trapped-ion results.
3. **7-site FMO:** Adolphs & Renger Hamiltonian. Main result — demonstrate scaling beyond 4-site limit.

**Metrics at each level:** Kraus rank vs χ, qubit count, circuit depth, trace distance from reference, wall-clock time (classical PT vs PT-VQD).

**Output:** Scaling tables and comparison plots (Makie.jl, consistent with classical-baselines styling). Data saved as structured files (CSV or JLD2) for reproducibility.

**Validation criteria:**
- Spin-boson and dimer match classical baselines within 1% trace distance
- Kraus rank stays bounded by χ as system size grows (core claim)
- 7-site FMO runs with ~11 qubits as predicted

**Common pitfalls:**
- FMO Hamiltonian has specific site energies and couplings from Adolphs 2006 — use exact values
- Each chromophore's bath has different reorganization energy
- Don't conflate "runs" with "is accurate" — always check convergence in Δt and χ

**Pivot signals** (from ideas report — flag these, don't just debug):
- Kraus rank is NOT bounded by χ in practice → undermines core scaling claim
- Stinespring circuit depth too large even for dimer → circuit approach impractical
- Trotter splitting introduces unacceptable errors at feasible Δt → fundamental method limitation

**Dependencies:** All previous skills (full pipeline).

**Trigger:** When running benchmark simulations, analyzing scaling results, or preparing comparison figures.

## Skill file locations

All skills will be placed in the repository's `.claude/skills/` directory:
- `.claude/skills/classical-baselines.md`
- `.claude/skills/pt-mpo-construction.md`
- `.claude/skills/kraus-extraction.md`
- `.claude/skills/quantum-circuit.md`
- `.claude/skills/benchmark-scaling.md`

## Design decisions

- **Lean over heavy:** Skills provide domain guardrails, not API tutorials. Claude can look up library docs; it can't look up physics constraints.
- **ITensors.jl for PT-MPO, not Yao.jl:** ITensors.jl has native MPO support with index management and SVD truncation. Output is converted to plain Julia arrays at the boundary so downstream stages don't depend on it.
- **Makie.jl for visualization:** Consistent publication-quality plots across all stages. PNG for inspection, PDF for paper.
- **Stage-specific triggers:** Each skill activates only when working on its pipeline stage, keeping context focused.
