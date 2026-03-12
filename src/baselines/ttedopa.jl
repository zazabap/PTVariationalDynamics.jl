using MPSDynamics

"""
    run_ttedopa(params::SpinBosonParams; dt=0.5, tfinal=20.0,
                chain_length=60, d=10, bond_dims=[10, 20, 40])

Run a T-TEDOPA simulation for the spin-boson model.

Returns a `SimulationResult` with observable "sigma_z".
Uses the highest bond dimension result from convergence comparison.

# MPSDynamics.jl API notes
- `chaincoeffs_finiteT` computes chain coefficients for finite-T TEDOPA
- `spinbosonmpo(ω0, Δ, d, N, cpars)` builds the MPO
  ω0 = epsilon (σ_z coefficient), Δ = tunneling (σ_x coefficient)
- `runsim` with `convobs` runs at multiple bond dimensions for convergence
- Data accessed via `dat["data/times"]`, `dat["convdata/sz"]`
"""
function run_ttedopa(params::SpinBosonParams;
                     dt::Float64=0.5,
                     tfinal::Float64=20.0,
                     chain_length::Int=60,
                     d::Int=10,
                     bond_dims::Vector{Int}=[10, 20, 40])
    N = chain_length

    # Chain coefficients for finite-temperature Ohmic bath (T-TEDOPA)
    # Third positional arg `ohmic=true` selects Ohmic spectral density
    # Keyword `α` (Greek) is the coupling strength (= ξ in our convention)
    cpars = chaincoeffs_finiteT(N, params.β, true;
        α=params.ξ, s=1, ωc=params.ωc,
        mc=4, mp=0, iq=1, idelta=2,
        procedure=:Lanczos, Mmax=5000, save=false
    )

    # Build MPO: H = (ε/2)σ_z + Δ σ_x + bath coupling
    # MPSDynamics spinbosonmpo takes (ω0, Δ, d, N, cpars)
    # where ω0 = gap (σ_z coeff), Δ = tunneling (σ_x coeff)
    H = spinbosonmpo(params.ε, params.Δ, d, N, cpars)

    # Initial state: spin-up ⊗ vacuum bath
    psi_up = unitcol(1, 2)
    A = productstatemps(physdims(H), state=[psi_up, fill(unitcol(1, d), N)...])

    # Observable: σ_z on site 1
    ob_sz = OneSiteObservable("sz", sz, 1)

    # Run with convergence comparison across bond dimensions
    A_final, dat = runsim(dt, tfinal, A, H;
        name="spin_boson_ttedopa",
        method=:TDVP1,
        obs=[],
        convobs=[ob_sz],
        convparams=bond_dims,
        verbose=false,
        save=false,
        plot=false
    )

    # Extract results — use highest bond dimension (last column)
    times = dat["data/times"]
    sz_data = dat["convdata/sz"]
    sigma_z = sz_data[:, end]  # last column = highest bond dim

    return SimulationResult(
        t=collect(times),
        observables=Dict("sigma_z" => collect(sigma_z)),
        params=params,
        method="T-TEDOPA",
        initial_state="up",
        metadata=Dict{String, Any}(
            "dt" => dt, "tfinal" => tfinal,
            "chain_length" => chain_length, "d" => d,
            "bond_dims" => bond_dims,
            "best_bond_dim" => bond_dims[end]
        )
    )
end
