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
