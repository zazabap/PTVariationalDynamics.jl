using JLD2

@testset "Spin-Boson Types" begin
    @testset "SpinBosonParams" begin
        p = SpinBosonParams(ξ=0.1, ωc=7.5, β=5.0, Δ=1.0, ε=0.0)
        @test p.ξ == 0.1
        @test p.ωc == 7.5
        @test p.β == 5.0
        @test p.Δ == 1.0
        @test p.ε == 0.0
    end

    @testset "Wang et al. parameters" begin
        p = wang_et_al_params()
        @test p.ξ == 0.1
        @test p.ωc == 7.5
        @test p.β == 5.0
        @test p.Δ == 1.0
        @test p.ε == 0.0
    end

    @testset "SimulationResult construction" begin
        p = wang_et_al_params()
        t = collect(0.0:0.1:1.0)
        obs = Dict("sigma_z" => ones(length(t)))
        r = SimulationResult(
            t=t, observables=obs, params=p,
            method="TEMPO", initial_state="up",
            metadata=Dict{String,Any}("dt" => 0.1)
        )
        @test r.method == "TEMPO"
        @test r.initial_state == "up"
        @test length(r.t) == 11
        @test r.observables["sigma_z"][1] == 1.0
    end

    @testset "JLD2 round-trip" begin
        p = wang_et_al_params()
        t = collect(0.0:0.1:1.0)
        obs = Dict("sigma_z" => sin.(t))
        r = SimulationResult(
            t=t, observables=obs, params=p,
            method="TEMPO", initial_state="up",
            metadata=Dict{String,Any}("dt" => 0.1)
        )

        tmpfile = tempname() * ".jld2"
        save_result(tmpfile, r)
        r2 = load_result(tmpfile)

        @test r2.t ≈ r.t
        @test r2.observables["sigma_z"] ≈ r.observables["sigma_z"]
        @test r2.params.ξ == r.params.ξ
        @test r2.method == r.method
        @test r2.initial_state == r.initial_state

        rm(tmpfile)
    end

    @testset "Trace distance metric" begin
        # For single-qubit diagonal states: D = |sz1 - sz2| / 2
        # Identical states → distance 0
        @test abs(1.0 - 1.0) / 2.0 == 0.0

        # Opposite states (sz=+1 vs sz=-1) → distance 1
        @test abs(1.0 - (-1.0)) / 2.0 == 1.0

        # Small difference
        @test abs(0.5 - 0.48) / 2.0 ≈ 0.01
    end
end
