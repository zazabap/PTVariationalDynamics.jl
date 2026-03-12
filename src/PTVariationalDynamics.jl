module PTVariationalDynamics

using LinearAlgebra
using JLD2
using CairoMakie

include("spectral_densities.jl")
include("spin_boson.jl")

# Baselines
include("baselines/tempo.jl")

# Plotting
include("plotting.jl")

export ohmic_spectral_density, drude_spectral_density
export SpinBosonParams, SimulationResult, wang_et_al_params
export save_result, load_result
export run_tempo
export plot_dynamics, plot_comparison, plot_convergence, save_figure

end
