using PTVariationalDynamics

base = wang_et_al_params()

# Sweep definitions: (name, parameter, values)
sweeps = [
    ("coupling_xi", :ξ, [0.01, 0.1, 0.5, 1.0]),
    ("temperature_beta", :β, [0.5, 1.0, 5.0, 20.0]),
    ("cutoff_wc", :ωc, [1.0, 5.0, 7.5, 20.0]),
]

mkpath("data/baselines")
mkpath("figures")

for (sweep_name, param_sym, values) in sweeps
    println("\n=== Sweep: $(sweep_name) (TEMPO only) ===")

    tempo_results = SimulationResult[]

    for val in values
        # Create modified params
        p = SpinBosonParams(;
            ξ = param_sym == :ξ ? val : base.ξ,
            ωc = param_sym == :ωc ? val : base.ωc,
            β = param_sym == :β ? val : base.β,
            Δ = base.Δ,
            ε = base.ε
        )

        println("  $(param_sym)=$(val): running TEMPO...")
        r_t = run_tempo(p; dt=0.1, ntimes=200, kmax=100)
        push!(tempo_results, r_t)
        save_result("data/baselines/spin_boson_tempo_$(sweep_name)_$(val).jld2", r_t)
    end

    # Plot TEMPO sweep
    fig = plot_convergence(tempo_results, string(param_sym), values)
    save_figure(fig, "figures/sweep_$(sweep_name)_tempo")

    println("  Figure saved to figures/sweep_$(sweep_name)_tempo.{png,pdf}")
end

println("\n=== All TEMPO sweeps complete ===")
println("\nNote: T-TEDOPA sweeps require separate environment.")
println("Run: julia --project=scripts/ttedopa_env scripts/parameter_sweep_ttedopa.jl")
