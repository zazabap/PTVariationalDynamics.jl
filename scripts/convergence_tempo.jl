using PTVariationalDynamics

params = wang_et_al_params()
mkpath("data/baselines")
mkpath("figures")

# --- Δt convergence ---
println("=== TEMPO Δt convergence ===")
dt_values = [0.1, 0.05, 0.025]
dt_results = SimulationResult[]

for dt in dt_values
    ntimes = round(Int, 20.0 / dt)  # t_final = 20
    println("  dt=$(dt), ntimes=$(ntimes)")
    r = run_tempo(params; dt=dt, ntimes=ntimes, kmax=100)
    push!(dt_results, r)
end

fig = plot_convergence(dt_results, "Δt", dt_values)
save_figure(fig, "figures/convergence_tempo_dt")

# Check: last two differ by < 0.1%
# Interpolate finer grid to coarser grid for comparison
using Interpolations
itp = linear_interpolation(dt_results[3].t, dt_results[3].observables["sigma_z"])
sz_fine_on_coarse = itp.(dt_results[2].t)
diff = maximum(abs.(dt_results[2].observables["sigma_z"] .- sz_fine_on_coarse))
println("Max diff between Δt=0.05 and Δt=0.025: $(diff)")

# --- K_max convergence ---
println("\n=== TEMPO K_max convergence ===")
kmax_values = [50, 100, 150, 200]
kmax_results = SimulationResult[]

for km in kmax_values
    println("  kmax=$(km)")
    r = run_tempo(params; dt=0.1, ntimes=200, kmax=km)
    push!(kmax_results, r)
end

fig = plot_convergence(kmax_results, "K_max", kmax_values)
save_figure(fig, "figures/convergence_tempo_kmax")

diff = maximum(abs.(kmax_results[3].observables["sigma_z"] .-
    kmax_results[4].observables["sigma_z"]))
println("Max diff between K_max=150 and K_max=200: $(diff)")

println("\nDone. Figures saved to figures/convergence_tempo_*.{png,pdf}")
