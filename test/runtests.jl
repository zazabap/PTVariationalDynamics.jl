using Test
using PTVariationalDynamics

@testset "PTVariationalDynamics" begin
    include("test_spectral_densities.jl")
    include("test_spin_boson.jl")
end
