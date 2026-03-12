# PT-VQD Claude Code Skills Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create five stage-specific Claude Code skills that guide implementation of the PTVariationalDynamics.jl computational pipeline.

**Architecture:** Each skill is a standalone markdown file in `.claude/skills/` with YAML frontmatter (name, description, trigger) and structured body (purpose, libraries, inputs/outputs, validation criteria, pitfalls). Skills reference each other by name for pipeline context but are independently invocable.

**Tech Stack:** Claude Code skills (markdown with YAML frontmatter), no code dependencies.

**Spec:** `docs/superpowers/specs/2026-03-12-ptvqd-skills-design.md`

---

## Chunk 1: Project setup and shared foundations

### Task 1: Create skills directory

**Files:**
- Create: `.claude/skills/` (directory)

- [ ] **Step 1: Create the skills directory**

Run: `mkdir -p .claude/skills`

- [ ] **Step 2: Commit**

```bash
git add -A .claude/skills
git commit -m "chore: create .claude/skills directory for PT-VQD pipeline skills"
```

### Task 2: Create `classical-baselines` skill

**Files:**
- Create: `.claude/skills/classical-baselines.md`

- [ ] **Step 1: Write the skill file**

```markdown
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
```

- [ ] **Step 2: Verify the file renders correctly**

Run: `head -5 .claude/skills/classical-baselines.md` to check frontmatter.

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/classical-baselines.md
git commit -m "feat: add classical-baselines Claude Code skill"
```

### Task 3: Create `pt-mpo-construction` skill

**Files:**
- Create: `.claude/skills/pt-mpo-construction.md`

- [ ] **Step 1: Write the skill file**

```markdown
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
```

- [ ] **Step 2: Verify frontmatter**

Run: `head -5 .claude/skills/pt-mpo-construction.md`

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/pt-mpo-construction.md
git commit -m "feat: add pt-mpo-construction Claude Code skill"
```

### Task 4: Create `kraus-extraction` skill

**Files:**
- Create: `.claude/skills/kraus-extraction.md`

- [ ] **Step 1: Write the skill file**

```markdown
---
name: kraus-extraction
description: Guide extraction of Kraus operators from PT-MPO slices via Choi matrix eigendecomposition with CPTP verification
---

# Kraus Extraction

## Purpose

Guide Claude when extracting Kraus operators from PT-MPO time-step slices via Choi matrix eigendecomposition. This is the bridge between the classical tensor network (PT-MPO) and the quantum circuit (Stinespring dilation).

## When to use

Invoke when writing Choi decomposition code, debugging Kraus operator extraction, or verifying CPTP (completely positive, trace-preserving) properties.

## Libraries

- `LinearAlgebra` (stdlib) — `eigen`, `Hermitian`, matrix operations. No heavy dependencies needed since inputs are small d²×d² matrices.
- `JLD2.jl` — serialization

## Inputs

A single time-step slice of the PT-MPO: a `d² × d²` matrix representing the CPTP map Φ_j for bath j. For a single qubit (d=2), this is a 4×4 matrix.

Loaded from `data/pt_mpo/{model}_bath{j}.jld2`, key `mpo_slices[n]`.

## Outputs

Save per-bath, per-step to `data/kraus/{model}_bath{j}_step{n}.jld2` with keys:
- `kraus_ops` — `Vector{Matrix{ComplexF64}}`, the Kraus operators
- `rank` — number of non-negligible Kraus operators
- `choi_eigenvalues` — full eigenvalue spectrum of the Choi matrix (for diagnostics)

## Algorithm

1. **Superoperator → Choi matrix:** Given the d²×d² superoperator S, form the Choi matrix:
   ```julia
   C = reshape(S, d, d, d, d)
   C = permutedims(C, (1, 3, 2, 4))
   C = reshape(C, d^2, d^2)
   ```
   Ensure C is Hermitian: `C = Hermitian((C + C') / 2)`

2. **Eigendecompose:** `λ, V = eigen(C)`. Eigenvalues should be real and ≥ 0.

3. **Extract Kraus operators:** For each eigenvalue λ_k > threshold:
   ```julia
   K_k = sqrt(λ_k) * reshape(V[:, k], d, d)
   ```

4. **Discard** eigenvalues below threshold (e.g., 1e-12).

## Validation criteria

- **CPTP completeness:** `sum(K' * K for K in kraus_ops) ≈ I` (use `norm(... - I) < 1e-10`)
- **Positive semidefinite Choi:** All eigenvalues ≥ 0 (after clamping noise)
- **Rank bound:** `rank ≤ χ` where χ is the PT-MPO bond dimension. This is the core claim of the paper.
- **Channel reconstruction:** `sum(K * ρ * K' for K in kraus_ops) ≈ S_reshaped * ρ_vec` for test input ρ

## Visualization

- Choi eigenvalue spectrum (bar plot) — shows rank structure
- Kraus rank vs χ across different bath parameters — verifies the rank bound

## Common pitfalls

- **Reshaping convention:** Julia is column-major. The `permutedims((1,3,2,4))` step is critical — getting this wrong produces a valid-looking but physically wrong Choi matrix. Always verify with the channel reconstruction test.
- **Negative eigenvalues from noise:** Numerical errors can produce small negative eigenvalues (e.g., -1e-15). Clamp these to zero, don't ignore them or treat them as a bug. But if eigenvalues are significantly negative (e.g., < -1e-6), the input superoperator is not CPTP — debug upstream.
- **Threshold choice:** Too aggressive → lose Kraus operators, channel is no longer trace-preserving. Too loose → include noise operators. Start with 1e-12 and verify CPTP completeness.

## Pipeline context

Stage 3 of 5. Depends on: `pt-mpo-construction` (for MPO slices). Feeds into: `quantum-circuit` (Kraus ops → Stinespring dilation).
```

- [ ] **Step 2: Verify frontmatter**

Run: `head -5 .claude/skills/kraus-extraction.md`

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/kraus-extraction.md
git commit -m "feat: add kraus-extraction Claude Code skill"
```

## Chunk 2: Circuit and benchmarking skills

### Task 5: Create `quantum-circuit` skill

**Files:**
- Create: `.claude/skills/quantum-circuit.md`

- [ ] **Step 1: Write the skill file**

```markdown
---
name: quantum-circuit
description: Guide Stinespring dilation circuit construction and p-VQD variational dynamics in Yao.jl for open quantum system simulation
---

# Quantum Circuit

## Purpose

Guide Claude when implementing Stinespring dilation circuits and variational quantum dynamics (p-VQD) in Yao.jl. This skill covers both converting Kraus operators into quantum circuits and running the variational time-stepping loop.

## When to use

Invoke when writing Yao.jl circuit code, implementing Stinespring dilation, debugging variational optimization, or during early Yao.jl warm-up/familiarization.

## Libraries

- `Yao.jl` — quantum circuit construction (`chain`, `put`, `control`, `matblock`), register management, measurement
- `YaoBlocks.jl` — custom gate blocks from matrices
- `Optim.jl` — BFGS fallback optimizer
- `JLD2.jl` — serialization
- `CairoMakie.jl` — result plots

## Inputs

- Kraus operators `{K_j^(k)}` loaded from `data/kraus/{model}_bath{j}_step{n}.jld2`
- System Hamiltonian parameters (coupling strengths, site energies)

## Outputs

Save to `data/circuit/{model}_dynamics.jld2` with keys:
- `t` — time array
- `density_matrices` — Vector of density matrices at each time step
- `expectation_values` — Dict of observable traces
- `circuit_params` — variational parameters θ at each step

## Sub-task 1: Stinespring dilation

Embed Kraus operators into a unitary on system + ancilla qubits.

**Construction:** Given r Kraus operators {K_1, ..., K_r} each d×d:
1. Ancilla size: `n_anc = ceil(Int, log2(r))` qubits. Pad r to next power of 2 if needed.
2. Build the isometry `V = [K_1; K_2; ...; K_r; zeros...]` of size `(d * 2^n_anc) × d`.
3. Extend V to a full unitary U (e.g., via QR completion or SVD).
4. Embed U as a gate in Yao.jl using `matblock(U)`.

**Two modes:**
- **(a) Mid-circuit measurement (default, ~11 qubits for FMO):** After applying U_j for bath j, measure ancilla qubits, reset to |0⟩, reuse for bath j+1. Use `measure!` + register reset.
- **(b) Deterministic dilation (~35 qubits for FMO):** Separate ancilla registers per bath, no measurement. Trace out all ancillas at the end. Use when mid-circuit measurement is unavailable or noisy.

Implement both. Default to (a).

## Sub-task 2: Variational ansatz (p-VQD)

Parameterized circuit V(θ_n) for system Hamiltonian evolution at each time step.

**Method:** Projected Variational Quantum Dynamics (Barison et al. 2021)
- **Cost function:** Fidelity F = |⟨ψ(θ_{n+1})|ψ_target⟩|² where ψ_target is the state after exact one-step evolution
- **Gradients:** Parameter-shift rule. For each parameter θ_i: ∂F/∂θ_i = [F(θ_i + π/2) - F(θ_i - π/2)] / 2
- **Convergence:** F > 1 - 1e-6
- **Fallback:** If p-VQD oscillates or fails to converge after 100 iterations, switch to direct fidelity maximization with BFGS via Optim.jl

## Time-stepping loop

At each time step t_n → t_{n+1}:
1. Apply variational circuit V(θ_n) for system Hamiltonian
2. For each bath j = 1, ..., N_sites: apply Stinespring circuit for Kraus operators at step n
3. Measure/trace-out ancillas
4. Record density matrix and observables
5. Optimize θ_{n+1} via p-VQD

## Validation criteria

- Stinespring output matches direct Kraus application: trace distance < 1e-8 (statevector sim)
- Ancilla in separable state after measurement/reset (no entanglement leakage)
- Full time-stepping loop reproduces `classical-baselines` dynamics for spin-boson and dimer
- Non-physical state check: after each step verify tr(ρ) = 1 and eigenvalues(ρ) ≥ 0. Apply trace renormalization and eigenvalue clamping if needed.

## Visualization

- Circuit fidelity F vs time step index
- Variational parameter θ convergence curves
- Ancilla measurement outcome statistics (should be dominated by |0⟩ for weak coupling)

## Common pitfalls

- **Qubit ordering:** Yao.jl uses LSB convention (qubit 1 is least significant bit). When mapping site indices to qubit indices, be explicit.
- **In-place measurement:** `measure!` collapses the register. Use `measure` (no bang) if you need the pre-measurement state. For density matrix simulation, use `DensityMatrix` register type.
- **Padding Kraus operators:** When r is not a power of 2, pad with zero matrices to fill the isometry. The unitary extension must still be valid (check unitarity: U†U = I).
- **Optimizer stalling:** Variational optimization can get stuck in local minima. Use multiple random initial θ (at least 3) and take the best. For early time steps, initialize θ near the identity circuit.

## Pipeline context

Stage 4 of 5. Depends on: `kraus-extraction` (Kraus operator inputs), `classical-baselines` (validation targets). Feeds into: `benchmark-scaling` (full pipeline results).
```

- [ ] **Step 2: Verify frontmatter**

Run: `head -5 .claude/skills/quantum-circuit.md`

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/quantum-circuit.md
git commit -m "feat: add quantum-circuit Claude Code skill"
```

### Task 6: Create `benchmark-scaling` skill

**Files:**
- Create: `.claude/skills/benchmark-scaling.md`

- [ ] **Step 1: Write the skill file**

```markdown
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
```

- [ ] **Step 2: Verify frontmatter**

Run: `head -5 .claude/skills/benchmark-scaling.md`

- [ ] **Step 3: Commit**

```bash
git add .claude/skills/benchmark-scaling.md
git commit -m "feat: add benchmark-scaling Claude Code skill"
```

### Task 7: Final verification and push

- [ ] **Step 1: Verify all 5 skill files exist**

Run: `ls -la .claude/skills/`

Expected: 5 `.md` files — classical-baselines, pt-mpo-construction, kraus-extraction, quantum-circuit, benchmark-scaling.

- [ ] **Step 2: Verify all frontmatter parses correctly**

Run: `head -4 .claude/skills/*.md`

Expected: Each file starts with `---`, has `name:` and `description:`, ends with `---`.

- [ ] **Step 3: Push to GitHub**

```bash
git push
```

- [ ] **Step 4: Verify on GitHub**

Run: `gh repo view zazabap/PTVariationalDynamics.jl --web` or check file listing with `gh api repos/zazabap/PTVariationalDynamics.jl/contents/.claude/skills`
