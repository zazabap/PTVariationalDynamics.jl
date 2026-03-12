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
