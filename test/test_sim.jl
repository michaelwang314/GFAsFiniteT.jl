using GFAsFiniteT

function run!()
    box = [5.0, 5.0, 5.0]
    num_steps = 10000000
    save_interval = trunc(Int64, num_steps / 10000)
    dt = 0.001
    kT = 1.0

    a = 1.0
    t, w = a / sqrt(2), 0.0

    north = [([0.0, a / 2, t / 2], "N1"), ([0.0, a / 2, -t / 2], "N2"), ([-w / 2, a / 2, 0.0], "N3"), ([w / 2, a / 2, 0.0], "N4")]
    south = [([0.0, -a / 2, t / 2], "S1"), ([0.0, -a / 2, -t / 2], "S2"), ([-w / 2, -a / 2, 0.0], "S3"), ([w / 2, -a / 2, 0.0], "S4")]
    east = [([a / 2, 0.0, t / 2], "E1"), ([a / 2, 0.0, -t / 2], "E2"), ([a / 2, w / 2, 0.0], "E3"), ([a / 2, -w / 2, 0.0], "E4")]
    west = [([-a / 2, 0.0, t / 2], "W1"), ([-a / 2, 0.0, -t / 2], "W2"), ([-a / 2, w / 2, 0.0], "W3"), ([-a / 2, -w / 2, 0.0], "W4")]
    linker_sites = [([a / 4, a / 4, 0.0], "linker_site"), ([-a / 4, a / 4, 0.0], "linker_site"), ([-a / 4, -a / 4, 0.0], "linker_site"), ([a / 4, -a / 4, 0.0], "linker_site")]
    excluders = [([0.0, 0.0, 0.0], "central")]
    body_particles = Vector{Particle}()
    for (r, id) in [north; south; east; west; linker_sites; excluders]
        push!(body_particles, Particle(r, id))
    end
    temp_rigid_body = RigidBody(body_particles, "")

    bodies = Vector{Body}()
    i = 0
    for Δr in [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [-0.5, -0.5, 0.0], [0.5, -0.5, 0.0]]
        new_rigid_body = deepcopy(temp_rigid_body)
        translate!(new_rigid_body, Δr .+ box ./ 2)
        set_body_ids!(new_rigid_body, "body $(i += 1)")
        push!(bodies, new_rigid_body)
    end
    
    all_particles = get_particle_list(bodies, [RigidBody])

    lj_particles = get_particles_with_ids(all_particles, ["central"])
    lj_r_cut = 2^(1 / 6) * a / sqrt(2)
    lj_imatrix = create_interaction_matrix([("central", "central"), ("central", "linker")])
    lj_cell_list = LinkedCellList(lj_particles, lj_r_cut, box)
    lj = LennardJones([("central", "central", Dict("ϵ" => 1.0, "σ" => a / sqrt(2), "r_cut" => lj_r_cut)), 
                       ("central", "linker", Dict("ϵ" => 1.0, "σ" => a / sqrt(2), "r_cut" => lj_r_cut))], lj_particles, lj_cell_list, lj_imatrix, box, false)

    attractor_imatrix = create_interaction_matrix([[("N$i", "S$i") for i = 1 : 4]; [("E$i", "W$i") for i = 1 : 4]])
    attractors = get_particles_with_ids(all_particles, ["$(d)$(i)" for d in ["N", "S", "E", "W"], i = 1 : 4][:])
    hb_bond_list = bind_closest(attractors, 0.01, attractor_imatrix)
    hb = HarmonicBond(25.0, 0.0, hb_bond_list, box, false)

    brownian = Brownian(bodies, dt, kT, box, false)

    system = System(bodies, [lj, hb], [lj_cell_list], brownian)
    trajectories = Trajectories(1, save_interval)
    run_simulation!(system, trajectories, num_steps)
    #save_system!(system, "TEST_OUTPUT/system.out")
    export_trajectories!(trajectories, "TEST_OUTPUT/trajectories.txt", [RigidBody])
end
run!()