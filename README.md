# PTVariationalDynamics.jl

Process-Tensor Variational Quantum Dynamics — a hybrid classical-quantum method for simulating non-Markovian multi-site open quantum systems.

## Overview

PT-VQD combines compressed process tensor MPOs (PT-MPOs) with variational quantum circuits to simulate open quantum system dynamics in regimes inaccessible to either approach alone:

- **Classical PT-MPO** efficiently captures non-Markovian bath memory, but hits a D⁴ outer bond bottleneck as the system grows.
- **Quantum circuits** represent multi-site system states on log₂(D) qubits, but need a way to incorporate bath effects.
- **PT-VQD** bridges the gap: Kraus operators extracted from compressed PT-MPOs (rank bounded by bond dimension χ) are implemented on the quantum circuit via Stinespring dilation.

```
Classical (PT-MPO)              Quantum Circuit
┌─────────────────────┐         ┌──────────────────────────┐
│ Per-bath PT-MPO     │         │ System qubits            │
│ via TEMPO / iTEBD   │─Kraus──▶│ + ⌈log₂χ⌉ ancilla qubits │
│ Bond dim: χ         │  ops    │ Variational ansatz V(θ)  │
└─────────────────────┘         └──────────────────────────┘
```

## Current Status

**Step 1: Classical Baselines** — in progress

| Component | Status |
|-----------|--------|
| TEMPO wrapper (QuantumDynamics.jl) | ✅ Done |
| T-TEDOPA wrapper (MPSDynamics.jl) | ✅ Done |
| Cross-validation (TEMPO vs T-TEDOPA) | ✅ Done |
| Convergence & parameter sweeps | 🔲 Scripts ready, not yet run |

**Steps 2–5** (PT-MPO extraction → Kraus pipeline → quantum circuit → FMO scaling) are planned but not yet started.

## Quick Start

### Prerequisites

- Julia 1.11+
- Git

### Installation

```bash
git clone https://github.com/zazabap/PTVariationalDynamics.jl.git
cd PTVariationalDynamics.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

### Run TEMPO Baseline

```bash
julia --project=. scripts/run_tempo_baseline.jl
```

### Run T-TEDOPA Baseline

T-TEDOPA uses a separate environment due to an ITensors version conflict:

```bash
julia --project=scripts/ttedopa_env -e 'using Pkg; Pkg.instantiate()'
julia --project=scripts/ttedopa_env scripts/run_ttedopa_baseline.jl
```

### Run Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Project Structure

```
src/
├── PTVariationalDynamics.jl   # Module root
├── spectral_densities.jl      # Ohmic and Drude spectral densities
├── spin_boson.jl              # SpinBosonParams, SimulationResult, serialization
├── plotting.jl                # CairoMakie plotting utilities
└── baselines/
    ├── tempo.jl               # TEMPO wrapper (QuantumDynamics.jl)
    └── ttedopa.jl             # T-TEDOPA wrapper (MPSDynamics.jl)

scripts/
├── run_tempo_baseline.jl      # Spin-boson TEMPO simulation
├── run_ttedopa_baseline.jl    # Spin-boson T-TEDOPA simulation
├── cross_validate.jl          # TEMPO vs T-TEDOPA comparison
├── convergence_tempo.jl       # TEMPO convergence (dt, kmax)
├── convergence_ttedopa.jl     # T-TEDOPA convergence (chain, bond dim)
├── parameter_sweep.jl         # TEMPO parameter sweep (ξ, β, ωc)
└── parameter_sweep_ttedopa.jl # T-TEDOPA parameter sweep

data/baselines/                # JLD2 simulation outputs
figures/                       # PNG + PDF figures
```

## Benchmark Parameters

Default spin-boson parameters follow [Wang et al. (ACS Omega, 2024)](https://doi.org/10.1021/acsomega.3c09720):

| Parameter | Symbol | Value |
|-----------|--------|-------|
| Coupling strength | ξ | 0.1 |
| Cutoff frequency | ωc | 7.5 |
| Inverse temperature | β | 5.0 |
| Tunneling | Δ | 1.0 |
| Bias | ε | 0.0 |

## Tech Stack

- **Tensor networks:** [ITensors.jl](https://github.com/ITensor/ITensors.jl)
- **Classical baselines:** [QuantumDynamics.jl](https://github.com/amartyabose/QuantumDynamics.jl) (TEMPO), [MPSDynamics.jl](https://github.com/shareloqs/MPSDynamics.jl) (T-TEDOPA)
- **Quantum circuits:** [Yao.jl](https://github.com/QuantumBFS/Yao.jl) (planned)
- **Visualization:** [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl)
- **Data format:** [JLD2](https://github.com/JuliaIO/JLD2.jl)

## License

MIT
