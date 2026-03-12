"""
    ohmic_spectral_density(ω; ξ, ωc)

Ohmic spectral density with exponential cutoff (Leggett convention):
J(ω) = 2π ξ ω exp(-ω/ωc)

where ξ is the dimensionless Kondo parameter and ωc is the cutoff frequency.
Note: some references use α = 2ξ or omit the 2π prefactor.
"""
function ohmic_spectral_density(ω; ξ, ωc)
    return 2π * ξ * ω * exp(-ω / ωc)
end

"""
    drude_spectral_density(ω; λ, γ)

Drude-Lorentz spectral density:
J(ω) = 2λ γ ω / (ω² + γ²)

where λ is the reorganization energy and γ is the characteristic frequency.
Used for the dimer and FMO models in later pipeline stages.
"""
function drude_spectral_density(ω; λ, γ)
    return 2λ * γ * ω / (ω^2 + γ^2)
end
