using PTVariationalDynamics

# Wang et al. 2024 parameters
params = wang_et_al_params()

println("Running TEMPO baseline...")
println("  ξ=$(params.ξ), ωc=$(params.ωc), β=$(params.β), Δ=$(params.Δ), ε=$(params.ε)")

result = run_tempo(params; dt=0.1, ntimes=200, kmax=100)

println("Done. Time range: [$(result.t[1]), $(result.t[end])]")
println("  ⟨σ_z(0)⟩ = $(result.observables["sigma_z"][1])")
println("  ⟨σ_z(end)⟩ = $(result.observables["sigma_z"][end])")

# Save
mkpath("data/baselines")
save_result("data/baselines/spin_boson_tempo.jld2", result)
println("Saved to data/baselines/spin_boson_tempo.jld2")
