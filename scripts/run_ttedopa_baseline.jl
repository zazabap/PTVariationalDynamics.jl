#!/usr/bin/env julia
# Run with: julia --project=scripts/ttedopa_env scripts/run_ttedopa_baseline.jl
#
# This script uses a separate Julia environment because MPSDynamics.jl
# requires ITensors 0.7.x, which conflicts with QuantumDynamics.jl (ITensors 0.9.x).
# Results are saved as JLD2 files compatible with the main project.

using MPSDynamics
using JLD2

# Wang et al. 2024 parameters
ξ = 0.1; ωc = 7.5; β = 5.0; Δ = 1.0; ε = 0.0

# T-TEDOPA parameters
dt = 0.5; tfinal = 20.0; chain_length = 60; d = 10
bond_dims = [10, 20, 40]

println("Running T-TEDOPA baseline...")
println("  ξ=$ξ, ωc=$ωc, β=$β, Δ=$Δ, ε=$ε")

N = chain_length

# Chain coefficients for finite-temperature Ohmic bath
cpars = chaincoeffs_finiteT(N, β;
    alpha=ξ, s=1, J=nothing, ωc=ωc,
    mc=4, mp=0, AB=nothing, iq=1, idelta=2,
    procedure=:Lanczos, Mmax=5000, save=false
)

# Build MPO
H = spinbosonmpo(ε, Δ, d, N, cpars)

# Initial state: spin-up ⊗ vacuum bath
psi_up = unitcol(1, 2)
A = productstatemps(physdims(H), state=[psi_up, fill(unitcol(1, d), N)...])

# Observable: σ_z on site 1
ob_sz = OneSiteObservable("sz", sz, 1)

# Run with convergence comparison
A_final, dat = runsim(dt, tfinal, A, H;
    name="spin_boson_ttedopa",
    method=:TDVP1,
    obs=[],
    convobs=[ob_sz],
    convparams=bond_dims,
    verbose=false,
    save=false,
    plot=false
)

# Extract results — use highest bond dimension (last column)
times = dat["data/times"]
sz_data = dat["convdata/sz"]
sigma_z = sz_data[:, end]

println("Done. Time range: [$(times[1]), $(times[end])]")
println("  ⟨σ_z(0)⟩ = $(sigma_z[1])")
println("  ⟨σ_z(end)⟩ = $(sigma_z[end])")

# Save in JLD2 format compatible with main project's load_result()
mkpath("data/baselines")
jldsave("data/baselines/spin_boson_ttedopa.jld2";
    t=collect(times),
    observables=Dict("sigma_z" => collect(sigma_z)),
    params=Dict{String,Any}("ξ"=>ξ, "ωc"=>ωc, "β"=>β, "Δ"=>Δ, "ε"=>ε),
    method="T-TEDOPA",
    initial_state="up",
    metadata=Dict{String, Any}(
        "dt" => dt, "tfinal" => tfinal,
        "chain_length" => chain_length, "d" => d,
        "bond_dims" => bond_dims,
        "best_bond_dim" => bond_dims[end]
    )
)
println("Saved to data/baselines/spin_boson_ttedopa.jld2")
