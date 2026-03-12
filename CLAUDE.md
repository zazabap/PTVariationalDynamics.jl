# PTVariationalDynamics.jl

Process Tensor Variational Quantum Dynamics — a hybrid classical-quantum method for simulating non-Markovian multi-site open quantum systems.

## Project overview

This Julia package implements the PT-VQD pipeline:
1. **Classical baselines** — reference dynamics via TEMPO and T-TEDOPA
2. **PT-MPO construction** — compress bath effects into process tensor MPOs (ITensors.jl)
3. **Kraus extraction** — Choi eigendecomposition of PT-MPO slices
4. **Quantum circuit** — Stinespring dilation and p-VQD in Yao.jl
5. **Benchmark scaling** — spin-boson → dimer → 7-site FMO

## Tech stack

- **Language:** Julia (1.11+)
- **Tensor networks:** ITensors.jl
- **Quantum circuits:** Yao.jl
- **Classical methods:** QuantumDynamics.jl, MPSDynamics.jl
- **Visualization:** CairoMakie.jl
- **Data format:** JLD2 (see `docs/superpowers/specs/2026-03-12-ptvqd-skills-design.md` for data conventions)

## Key files

- `articles/` — Typst ideas report and bibliography
- `src/` — Julia package source (future)
- `docs/superpowers/specs/` — design specs
- `docs/superpowers/plans/` — implementation plans
- `.claude/skills/` — 5 pipeline-stage skills for Claude Code

## Conventions

- Data outputs go to `data/{stage}/{model}_*.jld2`
- Figures go to `figures/` as PNG + PDF
- Compile Typst docs with `typst compile articles/<file>.typ`
- All PT-MPO outputs must be plain Julia arrays, not ITensor objects (boundary rule)
