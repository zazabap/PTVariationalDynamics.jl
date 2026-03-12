using Test
using PTVariationalDynamics

@testset "PTVariationalDynamics" begin
    include("test_spectral_densities.jl")
    include("test_spin_boson.jl")

    @testset "Integration smoke tests" begin
        if isdefined(PTVariationalDynamics, :run_tempo)
            @testset "TEMPO smoke test" begin
                p = wang_et_al_params()
                r = run_tempo(p; dt=0.5, ntimes=10, kmax=5)
                @test r.method == "TEMPO"
                @test length(r.t) > 0
                @test r.observables["sigma_z"][1] ≈ 1.0 atol=0.01
            end
        end
    end
end
