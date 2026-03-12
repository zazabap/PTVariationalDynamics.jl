"""
    SpinBosonParams(; ξ, ωc, β, Δ, ε)

Parameters for the spin-boson model.

Hamiltonian: H_s = (ε/2) σ_z + (Δ/2) σ_x
System-bath coupling: H_sb = σ_z ⊗ Σ_k g_k (a_k + a_k†)
Spectral density: J(ω) = 2π ξ ω exp(-ω/ωc)  (Ohmic, Leggett convention)
"""
Base.@kwdef struct SpinBosonParams
    ξ::Float64       # Kondo parameter (dimensionless coupling)
    ωc::Float64      # bath cutoff frequency
    β::Float64       # inverse temperature
    Δ::Float64       # tunneling (σ_x coefficient)
    ε::Float64       # bias (σ_z coefficient)
end

"""
    wang_et_al_params()

Return the benchmark parameters from Wang et al. 2024 (ACS Omega).
Unbiased spin-boson: ξ=0.1, ωc=7.5, β=5.0, Δ=1.0, ε=0.0.
"""
function wang_et_al_params()
    return SpinBosonParams(ξ=0.1, ωc=7.5, β=5.0, Δ=1.0, ε=0.0)
end

"""
    SimulationResult(; t, observables, params, method, initial_state, metadata)

Standardized output from any simulation method.
"""
Base.@kwdef struct SimulationResult
    t::Vector{Float64}
    observables::Dict{String, Vector{Float64}}
    params::SpinBosonParams
    method::String
    initial_state::String
    metadata::Dict{String, Any}
end

"""
    save_result(path, result::SimulationResult)

Save a SimulationResult to a JLD2 file.
"""
function save_result(path::String, result::SimulationResult)
    jldsave(path;
        t=result.t,
        observables=result.observables,
        params=result.params,
        method=result.method,
        initial_state=result.initial_state,
        metadata=result.metadata
    )
end

"""
    load_result(path) -> SimulationResult

Load a SimulationResult from a JLD2 file.
"""
function load_result(path::String)
    data = load(path)
    # Handle params saved as SpinBosonParams or as Dict (from separate T-TEDOPA env)
    p = data["params"]
    if p isa SpinBosonParams
        params = p
    else
        params = SpinBosonParams(
            ξ=p["ξ"], ωc=p["ωc"], β=p["β"], Δ=p["Δ"], ε=p["ε"]
        )
    end
    return SimulationResult(
        t=data["t"],
        observables=data["observables"],
        params=params,
        method=data["method"],
        initial_state=data["initial_state"],
        metadata=data["metadata"]
    )
end
