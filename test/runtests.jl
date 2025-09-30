using GFAsFiniteT
using Test

function test_bodies()
    try
        rigid_body = RigidBody([Particle([0.0, 0.0, 0.0], 1.0, "central"), 
                                Particle([1.0, 0.0, 0.0], 1.0, "east"), 
                                Particle([0.0, 1.0, 0.0], 1.0, "north"), 
                                Particle([-1.0, 0.0, 0.0], 1.0, "west"), 
                                Particle([0.0, -1.0, 0.0], 1.0, "south")], 
                                "body 1")
        println(rigid_body.centroid)
        println(rigid_body.γ_rot)
        println(rigid_body.axes)
        rotate!(rigid_body, [0.0, 0.0, 1.0], pi / 4)
        translate!(rigid_body, [10.0, 0.0, 0.0])
        println(rigid_body.particles)
        println(typeof(rigid_body.γ_rot))

        return true
    catch e
        println(e)
        return false
    end
end

function test_interactions()
    try
        rigid_bodies = Vector{RigidBody}()
        rigid_body = RigidBody([Particle([0.0, 0.0, 0.0], 1.0, "central"), 
                                Particle([0.5, 0.0, 0.0], 1.0, "east"), 
                                Particle([0.0, 0.5, 0.0], 1.0, "north"), 
                                Particle([-0.5, 0.0, 0.0], 1.0, "west"), 
                                Particle([0.0, -0.5, 0.0], 1.0, "south")], 
                                "")
        for Δr in [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [-0.5, -0.5, 0.0], [0.5, -0.5, 0.0]]
            new_rigid_body = deepcopy(rigid_body)
            translate!(new_rigid_body, Δr)

            push!(rigid_bodies, new_rigid_body)
        end


        return true
    catch e
        println(e)
        return false
    end
end

@testset "GFAsFiniteT.jl" begin
    #@test test_bodies()
    @test test_interactions()
end