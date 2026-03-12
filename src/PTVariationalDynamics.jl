module PTVariationalDynamics

using LinearAlgebra
using JLD2

include("spectral_densities.jl")
include("spin_boson.jl")

export ohmic_spectral_density, drude_spectral_density
export SpinBosonParams, SimulationResult, wang_et_al_params
export save_result, load_result

end
