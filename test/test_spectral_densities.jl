@testset "Spectral Densities" begin
    @testset "Ohmic" begin
        # J(ω) = 2π ξ ω exp(-ω/ωc)  (Leggett convention)
        ξ = 0.1; ωc = 7.5

        # At ω = 0, J = 0
        @test ohmic_spectral_density(0.0; ξ=ξ, ωc=ωc) == 0.0

        # At ω = ωc, J = 2π ξ ωc exp(-1)
        expected = 2π * ξ * ωc * exp(-1)
        @test ohmic_spectral_density(ωc; ξ=ξ, ωc=ωc) ≈ expected

        # Positive for ω > 0
        @test ohmic_spectral_density(1.0; ξ=ξ, ωc=ωc) > 0.0

        # Linearity in ξ
        J1 = ohmic_spectral_density(1.0; ξ=0.1, ωc=ωc)
        J2 = ohmic_spectral_density(1.0; ξ=0.2, ωc=ωc)
        @test J2 ≈ 2 * J1
    end

    @testset "Drude" begin
        # J(ω) = 2λ γ ω / (ω² + γ²)
        λ = 35.0; γ = 106.18

        # At ω = 0, J = 0
        @test drude_spectral_density(0.0; λ=λ, γ=γ) == 0.0

        # At ω = γ, J = 2λ γ² / (2γ²) = λ
        @test drude_spectral_density(γ; λ=λ, γ=γ) ≈ λ

        # Positive for ω > 0
        @test drude_spectral_density(1.0; λ=λ, γ=γ) > 0.0
    end
end
