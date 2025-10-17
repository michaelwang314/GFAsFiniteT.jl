using GFAsFiniteT

function run!()
    box = [5.0, 5.0, 5.0]
    num_steps = 10000000
    save_interval = trunc(Int64, num_steps / 10000)
    dt = 0.0001
    kT = 0.1

    γ_linker = 1.0;
    γ_constituents = 1.0;

    a = 1.0
    t, w = a / sqrt(2), 0.0
    k_corners = 10.0

    num_linkers = 10

    north = [([0.0, a / 2, t / 2], γ_constituents, "N1"), ([0.0, a / 2, -t / 2], γ_constituents, "N2"), ([-w / 2, a / 2, 0.0], γ_constituents, "N3"), ([w / 2, a / 2, 0.0], γ_constituents, "N4")]
    south = [([0.0, -a / 2, t / 2], γ_constituents, "S1"), ([0.0, -a / 2, -t / 2], γ_constituents, "S2"), ([-w / 2, -a / 2, 0.0], γ_constituents, "S3"), ([w / 2, -a / 2, 0.0], γ_constituents, "S4")]
    east = [([a / 2, 0.0, t / 2], γ_constituents, "E1"), ([a / 2, 0.0, -t / 2], γ_constituents, "E2"), ([a / 2, w / 2, 0.0], γ_constituents, "E3"), ([a / 2, -w / 2, 0.0], γ_constituents, "E4")]
    west = [([-a / 2, 0.0, t / 2], γ_constituents, "W1"), ([-a / 2, 0.0, -t / 2], γ_constituents, "W2"), ([-a / 2, w / 2, 0.0], γ_constituents, "W3"), ([-a / 2, -w / 2, 0.0], γ_constituents, "W4")]
    linker_sites = [([a / 4, a / 4, 0.0], γ_constituents, "linker_site"), ([-a / 4, a / 4, 0.0], γ_constituents, "linker_site"), ([-a / 4, -a / 4, 0.0], γ_constituents, "linker_site"), ([a / 4, -a / 4, 0.0], γ_constituents, "linker_site")]
    excluders = [([0.0, 0.0, 0.0], γ_constituents, "central")]
    body_particles = Vector{Particle}()
    for (r, γ, id) in [north; south; east; west; linker_sites; excluders]
        push!(body_particles, Particle(r, γ, id))
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

    for n = 1 : num_linkers
        push!(bodies, Particle(5 * rand(3), γ_linker, "linker"))
    end
    
    all_particles = get_particle_list(bodies)

    r_linker_site = 0.125
    r_linker = 0.125
    m_r_cut = r_linker_site + r_linker
    m_particles = get_particles_with_ids(all_particles, ["linker_site", "linker"])
    m_cell_list = LinkedCellList(m_particles, m_r_cut, box)
    m_params = [("linker_site", "linker", Dict("D0" => 1.0, "α" => 2 / m_r_cut, "r0" => 0.0, "r_cut" => m_r_cut))]
    m = Morse(m_params, m_particles, m_cell_list, box; multithreaded = false)

    r_excluder = a / (2 * sqrt(2))
    lj_r_cut = 2^(1 / 6) * 2 * r_excluder
    lj_particles = get_particles_with_ids(all_particles, ["central", "linker"])
    lj_cell_list = LinkedCellList(lj_particles, lj_r_cut, box)
    lj_params = [("central", "central", Dict("ϵ" => 1.0, "σ" => 2 * r_excluder, "r_cut" => lj_r_cut)),
                 ("linker", "linker", Dict("ϵ" => 1.0, "σ" => 2 * r_linker, "r_cut" => 2^(1 / 6) * 2 * r_linker))]
    lj = LennardJones(lj_params, lj_particles, lj_cell_list, box; multithreaded = false)

    attractor_imatrix = create_interaction_matrix([[("N$i", "S$i") for i = 1 : 4]; [("E$i", "W$i") for i = 1 : 4]])
    attractors = get_particles_with_ids(all_particles, ["$(d)$(i)" for d in ["N", "S", "E", "W"], i = 1 : 4][:])
    hb_bond_list = bind_closest(attractors, 0.01, attractor_imatrix)
    hb = HarmonicBond(k_corners, 0.0, hb_bond_list, box; multithreaded = false)

    brownian = Brownian(bodies, dt, kT, box; multithreaded = false)

    system = System(bodies, [lj, hb, m], [lj_cell_list, m_cell_list], brownian)
    trajectories = Trajectories(save_interval)
    run_simulation!(system, trajectories, num_steps)
    #save_system!(system, "TEST_OUTPUT/system.out")
    export_trajectories!(trajectories, "TEST_OUTPUT/trajectories_tetramer_test.txt")
end
run!()