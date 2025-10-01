using GFAsFiniteT
using Test

function test_bodies()
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
end

function test_cell_lists()
    rigid_bodies = Vector{RigidBody}()
    rigid_body = RigidBody([Particle([0.0, 0.0, 0.0], 1.0, "central"), 
                            Particle([0.5, 0.0, 0.0], 1.0, "east"), 
                            Particle([0.0, 0.5, 0.0], 1.0, "north"), 
                            Particle([-0.5, 0.0, 0.0], 1.0, "west"), 
                            Particle([0.0, -0.5, 0.0], 1.0, "south")], 
                            "")
    body_id = 1
    for Δr in [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [-0.5, -0.5, 0.0], [0.5, -0.5, 0.0]]
        new_rigid_body = deepcopy(rigid_body)
        set_rigid_body_id!(new_rigid_body, "body $(body_id)")
        translate!(new_rigid_body, Δr .+ [2.0, 2.0, 2.0])

        push!(rigid_bodies, new_rigid_body)
        body_id += 1
    end
    particles = rigid_bodies_to_particle_list(rigid_bodies)

    cell_list = LinkedCellList(particles, [1.0, 1.0, 1.0], [4.0, 4.0, 4.0])
    update_cell_list!(cell_list)
    bond_list = bind_closest(particles, 0.01, create_interaction_matrix([("north", "south"), ("east", "west")]))
    for bond in bond_list.bonds
        println(bond)
    end
end

@testset "GFAsFiniteT.jl" begin
    #test_bodies()
    test_cell_lists()
end