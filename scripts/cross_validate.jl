using PTVariationalDynamics
using LinearAlgebra
using Interpolations

"""
    trace_distance_from_sigma_z(sz1, sz2)

Compute trace distance between two single-qubit states from their σ_z expectation values.
For a single qubit, ρ = (I + ⟨σ⟩·σ)/2, and the trace distance simplifies when only σ_z differs.
D(ρ₁, ρ₂) = ½ tr|ρ₁ - ρ₂|. For diagonal states (σ_z only): D = |sz1 - sz2| / 2.
"""
function trace_distance_from_sigma_z(sz1::Float64, sz2::Float64)
    return abs(sz1 - sz2) / 2.0
end

println("=== Cross-Validation: TEMPO vs T-TEDOPA ===")

# Run TEMPO (or load if already saved)
params = wang_et_al_params()
println("Parameters: ξ=$(params.ξ), ωc=$(params.ωc), β=$(params.β)")

tempo_path = "data/baselines/spin_boson_tempo.jld2"
ttedopa_path = "data/baselines/spin_boson_ttedopa.jld2"

if isfile(tempo_path)
    println("\nLoading TEMPO from $tempo_path")
    r_tempo = load_result(tempo_path)
else
    println("\nRunning TEMPO...")
    r_tempo = run_tempo(params; dt=0.1, ntimes=200, kmax=100)
    mkpath("data/baselines")
    save_result(tempo_path, r_tempo)
end

if isfile(ttedopa_path)
    println("Loading T-TEDOPA from $ttedopa_path")
    r_ttedopa = load_result(ttedopa_path)
else
    error("T-TEDOPA data not found at $ttedopa_path. Run first:\n  julia --project=scripts/ttedopa_env scripts/run_ttedopa_baseline.jl")
end

# Interpolate T-TEDOPA to TEMPO time grid for comparison
t_tempo = r_tempo.t
sz_tempo = r_tempo.observables["sigma_z"]

itp = linear_interpolation(r_ttedopa.t, r_ttedopa.observables["sigma_z"])
sz_ttedopa_interp = itp.(t_tempo)

# Compute metrics
pointwise_diff = abs.(sz_tempo .- sz_ttedopa_interp)
max_diff = maximum(pointwise_diff)
trace_dists = [trace_distance_from_sigma_z(sz_tempo[i], sz_ttedopa_interp[i]) for i in eachindex(t_tempo)]
max_trace_dist = maximum(trace_dists)

println("\n=== Results ===")
println("Max |Δ⟨σ_z⟩|: $(round(max_diff; digits=6))")
println("Max trace distance: $(round(max_trace_dist; digits=6))")

if max_trace_dist < 0.01
    println("✓ PASSED: Trace distance < 0.01 (1% threshold)")
else
    println("✗ FAILED: Trace distance $(max_trace_dist) ≥ 0.01")
    println("  Debug before proceeding to parameter sweeps.")
end

# Plot
mkpath("figures")
fig = plot_comparison(r_tempo,
    SimulationResult(
        t=t_tempo,
        observables=Dict("sigma_z" => sz_ttedopa_interp),
        params=params, method="T-TEDOPA (interp)",
        initial_state="up", metadata=Dict{String,Any}()
    )
)
save_figure(fig, "figures/cross_validation_sigma_z")
println("\nFigure saved to figures/cross_validation_sigma_z.{png,pdf}")
