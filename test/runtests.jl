using GFAsFiniteT
using Test

@testset "GFAsFiniteT.jl" begin
    include("test_sim.jl")
    include("interaction_test.jl")
end