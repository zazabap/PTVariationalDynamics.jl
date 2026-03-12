using QuantumDynamics

"""
    run_tempo(params::SpinBosonParams; dt=0.1, ntimes=200, kmax=100,
              svec=[1.0 -1.0], cutoff=1e-13, maxdim=500)

Run a TEMPO simulation for the spin-boson model.

Returns a `SimulationResult` with observables "sigma_z" and "sigma_x".

# QuantumDynamics.jl API notes
- `Utilities.create_tls_hamiltonian(; ϵ, Δ)` builds H = (ε/2)σ_z - (Δ/2)σ_x
  Note the minus sign on Δ: our convention H_s = (ε/2)σ_z + (Δ/2)σ_x requires Δ_QD = -Δ
- `SpectralDensities.ExponentialCutoff(; ξ, ωc)` matches our Leggett convention
- `svec = [1.0 -1.0]` encodes σ_z coupling: diagonal elements ⟨↑|σ_z|↑⟩=1, ⟨↓|σ_z|↓⟩=-1
"""
function run_tempo(params::SpinBosonParams;
                   dt::Float64=0.1,
                   ntimes::Int=200,
                   kmax::Int=100,
                   svec=[1.0 -1.0],
                   cutoff::Float64=1e-13,
                   maxdim::Int=500)
    # Build system Hamiltonian
    # QuantumDynamics uses H = (ε/2)σ_z - (Δ/2)σ_x, so pass -Δ to get +Δ σ_x
    H = Utilities.create_tls_hamiltonian(; ϵ=params.ε, Δ=-params.Δ)

    # Spectral density (Ohmic with exponential cutoff)
    Jw = SpectralDensities.ExponentialCutoff(; ξ=params.ξ, ωc=params.ωc)

    # Initial state: |↑⟩ = [1, 0] → ρ₀ = |↑⟩⟨↑|
    ρ0 = ComplexF64[1.0 0.0; 0.0 0.0]

    # Bare forward-backward propagators
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)

    # Run TEMPO
    extraargs = TEMPO.TEMPOArgs(; cutoff=cutoff, maxdim=maxdim)
    t, ρs = TEMPO.propagate(;
        fbU=fbU, Jw=[Jw], β=params.β, ρ0=ρ0,
        dt=dt, ntimes=ntimes, kmax=kmax,
        extraargs=extraargs, svec=svec
    )

    # Extract observables from density matrices
    # ρs[:, :, k] is the 2×2 density matrix at time t[k]
    nsteps = length(t)
    sigma_z = Vector{Float64}(undef, nsteps)
    sigma_x = Vector{Float64}(undef, nsteps)
    for k in 1:nsteps
        sigma_z[k] = real(ρs[1, 1, k] - ρs[2, 2, k])
        sigma_x[k] = real(ρs[1, 2, k] + ρs[2, 1, k])
    end

    return SimulationResult(
        t=collect(t),
        observables=Dict("sigma_z" => sigma_z, "sigma_x" => sigma_x),
        params=params,
        method="TEMPO",
        initial_state="up",
        metadata=Dict{String, Any}(
            "dt" => dt, "ntimes" => ntimes, "kmax" => kmax,
            "cutoff" => cutoff, "maxdim" => maxdim
        )
    )
end
