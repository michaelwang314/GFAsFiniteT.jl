using GFAsFiniteT

function run!()
    box = [5.0, 5.0, 5.0]
    num_steps = 1000000
    save_interval = trunc(Int64, num_steps / 10000)
    dt = 0.0001
    kT = 1.0

    a = 1.0
    t, w = a / sqrt(2), 0.0
    k_corners = 1000.0

    num_linkers = 0

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
    temp_rigid_body = RigidBody(body_particles)

    bodies = Vector{Body}()
    i = 0
    for Δr in [[0.5, 0.5, 0.0], [-0.5, 0.5, 0.0], [-0.5, -0.5, 0.0], [0.5, -0.5, 0.0]]
        new_rigid_body = deepcopy(temp_rigid_body)
        translate!(new_rigid_body, Δr .+ box ./ 2)
        set_body_ids!(new_rigid_body, "body $(i += 1)")
        push!(bodies, new_rigid_body)
    end

    # add linkers
    push!(bodies, Particle(box ./ 2, "linker"))
    push!(bodies, Particle([1.0, 1.0, 1.0], "linker"))

    all_particles = get_particle_list(bodies, [RigidBody])

    r_excluder = a / (2 * sqrt(2))
    r_linker = 0.1
    lj_r_cut = maximum(2^(1 / 6) * [2 * r_excluder, r_excluder + r_linker, 2 * r_linker])
    lj_particles = get_particles_with_ids(all_particles, ["central", "linker"])
    lj_cell_list = LinkedCellList(lj_particles, lj_r_cut, box)
    lj_params = [("central", "central", Dict("ϵ" => 1.0, "σ" => 2 * r_excluder, "r_cut" => 2^(1 / 6) * 2 * r_excluder)), 
                 ("central", "linker", Dict("ϵ" => 1.0, "σ" => r_excluder + r_linker, "r_cut" => 2^(1 / 6) * (r_excluder + r_linker))),
                 ("linker", "linker", Dict("ϵ" => 1.0, "σ" => 2 * r_linker, "r_cut" => 2^(1 / 6) * 2 * r_linker))]
    lj = LennardJones(lj_params, lj_particles, lj_cell_list, box, false)

    attractor_imatrix = create_interaction_matrix([[("N$i", "S$i") for i = 1 : 4]; [("E$i", "W$i") for i = 1 : 4]])
    attractors = get_particles_with_ids(all_particles, ["$(d)$(i)" for d in ["N", "S", "E", "W"], i = 1 : 4][:])
    hb_bond_list = bind_closest(attractors, 0.01, attractor_imatrix)
    hb = HarmonicBond(k_corners, 0.0, hb_bond_list, box, false)

    r_linker_site = 0.1
    m_r_cut = r_linker + r_linker_site
    m_particles = get_particles_with_ids(all_particles, ["linker", "linker_site"])
    m_cell_list = LinkedCellList(m_particles, m_r_cut, box)
    m_params = [("linker", "linker_site", Dict("D0" => 15.0, "α" => 50.0, "r0" => 0.0, "r_cut" => m_r_cut))]
    m = Morse(m_params, m_particles, m_cell_list, box, false)

    brownian = Brownian(bodies, dt, kT, box, false)

    system = System(bodies, [lj, m, hb], [lj_cell_list], brownian)
    trajectories = Trajectories(1, save_interval)
    run_simulation!(system, trajectories, num_steps)
    #save_system!(system, "TEST_OUTPUT/system.out")
    export_trajectories!(trajectories, "TEST_OUTPUT/trajectories.txt", [RigidBody])
end
run!()