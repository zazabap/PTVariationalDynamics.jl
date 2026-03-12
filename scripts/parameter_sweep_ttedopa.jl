#!/usr/bin/env julia
# Run with: julia --project=scripts/ttedopa_env scripts/parameter_sweep_ttedopa.jl
#
# Uses separate environment due to ITensors version conflict.

using MPSDynamics
using JLD2

function run_ttedopa_standalone(; ξ, ωc, β, Δ, ε, dt=0.5, tfinal=20.0,
                                  chain_length=60, d=10, bond_dims=[10, 20, 40])
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
        name="ttedopa_sweep", method=:TDVP1,
        obs=[], convobs=[ob_sz], convparams=bond_dims,
        verbose=false, save=false, plot=false
    )
    times = dat["data/times"]
    sz_data = dat["convdata/sz"]
    return collect(times), collect(sz_data[:, end])
end

# Base parameters
base_ξ = 0.1; base_ωc = 7.5; base_β = 5.0; base_Δ = 1.0; base_ε = 0.0

sweeps = [
    ("coupling_xi", "ξ", [(0.01, base_ωc, base_β), (0.1, base_ωc, base_β),
                           (0.5, base_ωc, base_β), (1.0, base_ωc, base_β)]),
    ("temperature_beta", "β", [(base_ξ, base_ωc, 0.5), (base_ξ, base_ωc, 1.0),
                                (base_ξ, base_ωc, 5.0), (base_ξ, base_ωc, 20.0)]),
    ("cutoff_wc", "ωc", [(base_ξ, 1.0, base_β), (base_ξ, 5.0, base_β),
                           (base_ξ, 7.5, base_β), (base_ξ, 20.0, base_β)]),
]

mkpath("data/baselines")

for (sweep_name, param_name, param_sets) in sweeps
    println("\n=== Sweep: $(sweep_name) (T-TEDOPA) ===")

    for (ξ_val, ωc_val, β_val) in param_sets
        val = param_name == "ξ" ? ξ_val : param_name == "ωc" ? ωc_val : β_val
        println("  $(param_name)=$(val): running T-TEDOPA...")
        t, sigma_z = run_ttedopa_standalone(; ξ=ξ_val, ωc=ωc_val, β=β_val,
                                              Δ=base_Δ, ε=base_ε)
        jldsave("data/baselines/spin_boson_ttedopa_$(sweep_name)_$(val).jld2";
            t=t, observables=Dict("sigma_z" => sigma_z),
            params=Dict{String,Any}("ξ"=>ξ_val, "ωc"=>ωc_val, "β"=>β_val,
                                     "Δ"=>base_Δ, "ε"=>base_ε),
            method="T-TEDOPA", initial_state="up",
            metadata=Dict{String,Any}("sweep"=>sweep_name, "param_value"=>val)
        )
    end
end

println("\n=== All T-TEDOPA sweeps complete ===")
println("Data saved to data/baselines/spin_boson_ttedopa_*.jld2")
println("Run plotting from main project to generate figures.")
