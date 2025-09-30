using GFAsFiniteT
using Test

@testset "GFAsFiniteT.jl" begin
    particles = [Particle([0.0, 0.0, 0.0], 1.0, 5, 1), Particle([1.0, 0.0, 0.0], 1.0, 1, 1), Particle([0.0, 1.0, 0.0], 1.0, 2, 1), Particle([-1.0, 0.0, 0.0], 1.0, 3, 1), Particle([0.0, -1.0, 0.0], 1.0, 4, 1)]
    rigid_body = RigidBody(particles, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [1.0, 1.0, 1.0], 1.0)
    println(rigid_body.centroid)
    println(rigid_body.γ_rot)
    println(rigid_body.axes)
    rotate!(rigid_body, [0.0, 0.0, 1.0], pi / 4)
    println(rigid_body.particles)
    println(typeof(rigid_body.γ_rot))
end
