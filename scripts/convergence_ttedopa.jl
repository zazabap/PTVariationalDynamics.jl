#!/usr/bin/env julia
# Run with: julia --project=scripts/ttedopa_env scripts/convergence_ttedopa.jl
#
# Uses separate environment due to ITensors version conflict.
# Saves JLD2 files loadable by main project.

using MPSDynamics
using JLD2

# Wang et al. 2024 parameters
ξ = 0.1; ωc = 7.5; β = 5.0; Δ = 1.0; ε = 0.0

function run_ttedopa_standalone(; ξ, ωc, β, Δ, ε, dt=0.5, tfinal=20.0,
                                  chain_length=60, d=10, bond_dims=[40])
    N = chain_length
    cpars = chaincoeffs_finiteT(N, β, true;
        α=ξ, s=1, ωc=ωc,
        mc=4, mp=0, iq=1, idelta=2,
        procedure=:Lanczos, Mmax=500, save=false
    )
    H = spinbosonmpo(ε, Δ, d, N, cpars)
    psi_up = unitcol(1, 2)
    A = productstatemps(physdims(H), state=[psi_up, fill(unitcol(1, d), N)...])
    ob_sz = OneSiteObservable("sz", sz, 1)
    A_final, dat = runsim(dt, tfinal, A, H;
        name="ttedopa_conv", method=:TDVP1,
        obs=[], convobs=[ob_sz], convparams=bond_dims,
        verbose=false, save=false, plot=false
    )
    times = dat["data/times"]
    sz_data = dat["convdata/sz"]
    return collect(times), collect(sz_data[:, end])
end

mkpath("data/baselines")

# --- Chain length convergence ---
println("=== T-TEDOPA chain length convergence ===")
chain_values = [30, 60, 90, 120]

for cl in chain_values
    println("  chain_length=$cl")
    t, sigma_z = run_ttedopa_standalone(; ξ=ξ, ωc=ωc, β=β, Δ=Δ, ε=ε,
                                          chain_length=cl, bond_dims=[40])
    jldsave("data/baselines/conv_ttedopa_chain_$cl.jld2";
        t=t, observables=Dict("sigma_z" => sigma_z),
        params=Dict{String,Any}("ξ"=>ξ, "ωc"=>ωc, "β"=>β, "Δ"=>Δ, "ε"=>ε),
        method="T-TEDOPA", initial_state="up",
        metadata=Dict{String,Any}("chain_length"=>cl, "bond_dim"=>40)
    )
end

if length(chain_values) >= 2
    # Quick convergence check for last two
    d1 = load("data/baselines/conv_ttedopa_chain_$(chain_values[end-1]).jld2")
    d2 = load("data/baselines/conv_ttedopa_chain_$(chain_values[end]).jld2")
    diff = maximum(abs.(d1["observables"]["sigma_z"] .- d2["observables"]["sigma_z"]))
    println("Max diff between N=$(chain_values[end-1]) and N=$(chain_values[end]): $diff")
end

# --- Bond dimension convergence ---
println("\n=== T-TEDOPA bond dimension convergence ===")
bdim_values = [10, 20, 40, 60]

for bd in bdim_values
    println("  bond_dim=$bd")
    t, sigma_z = run_ttedopa_standalone(; ξ=ξ, ωc=ωc, β=β, Δ=Δ, ε=ε,
                                          chain_length=60, bond_dims=[bd])
    jldsave("data/baselines/conv_ttedopa_bdim_$bd.jld2";
        t=t, observables=Dict("sigma_z" => sigma_z),
        params=Dict{String,Any}("ξ"=>ξ, "ωc"=>ωc, "β"=>β, "Δ"=>Δ, "ε"=>ε),
        method="T-TEDOPA", initial_state="up",
        metadata=Dict{String,Any}("chain_length"=>60, "bond_dim"=>bd)
    )
end

if length(bdim_values) >= 2
    d1 = load("data/baselines/conv_ttedopa_bdim_$(bdim_values[end-1]).jld2")
    d2 = load("data/baselines/conv_ttedopa_bdim_$(bdim_values[end]).jld2")
    diff = maximum(abs.(d1["observables"]["sigma_z"] .- d2["observables"]["sigma_z"]))
    println("Max diff between D=$(bdim_values[end-1]) and D=$(bdim_values[end]): $diff")
end

println("\nDone. Data saved to data/baselines/conv_ttedopa_*.jld2")
println("Run plotting from main project to generate figures.")
