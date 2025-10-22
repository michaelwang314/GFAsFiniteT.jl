using GFAsFiniteT

function run!()
    box = [20.0, 20.0, 20.0]
    num_steps = 1000000
    save_interval = trunc(Int64, num_steps / 10000)
    dt = 0.0001
    kT = 0.1

    γ = 1.0

    # Pair of particles interacting via LJ
    ϵ = 1.0
    σ = 1.0
    lj_r_cut = 2 * 2^(1 / 6) * σ
    lj_params = [(:lj_particle, :lj_particle, Dict(:ϵ => ϵ, :σ => σ, :r_cut => lj_r_cut))]
    lj_particle_pair = [Particle(box ./ 2 .- [lj_r_cut / 4, 0.0, 0.0], γ, :lj_particle), Particle(box ./ 2 .+ [lj_r_cut / 4, 0.0, 0.0], γ, :lj_particle)]
    lj_cell_list, _ = LinkedCellList(lj_particle_pair, lj_r_cut, box)
    lj = LennardJones(lj_params, lj_cell_list, box; multithreaded = false)

    # Pair of particles interacting via Morse
    D0 = 5.0
    α = 10.0
    r0 = 0.1
    m_r_cut = 2.0
    m_params = [(:morse_particle, :morse_particle, Dict(:D0 => D0, :α => α, :r0 => r0, :r_cut => m_r_cut))]
    m_particle_pair = [Particle(box ./ 2, γ, :morse_particle), Particle(box ./ 2, γ, :morse_particle)]
    m_cell_list, _ = LinkedCellList(m_particle_pair, m_r_cut, box)
    m = Morse(m_params, m_cell_list, box; multithreaded = false)

    # Pair of particles interacting via a harmonic bond
    k = 20.0
    r0 = 0.5
    hb_particle_pair = [Particle(box ./ 2 .- [r0 / 2, 0.0, 0.0], γ, :harmonic_particle), Particle(box ./ 2 .+ [r0 / 2, 0.0, 0.0], γ, :harmonic_particle)]
    hb_imatrix = create_interaction_matrix([(:harmonic_particle, :harmonic_particle)])
    hb_bond_list = bind_closest(hb_particle_pair, 2 * r0, hb_imatrix)
    hb = HarmonicBond(k, r0, hb_bond_list, box; multithreaded = false)

    # Pair of rigid bodies interacting
    k_corners = 500.0
    a, t, w = 1.0, 1.0 / sqrt(2), 0.25 / sqrt(2)
    north = [([0.0, a / 2, t / 2], γ, :N1), ([0.0, a / 2, -t / 2], γ, :N2), ([-w / 2, a / 2, 0.0], γ, :N3), ([w / 2, a / 2, 0.0], γ, :N4)]
    south = [([0.0, -a / 2, t / 2], γ, :S1), ([0.0, -a / 2, -t / 2], γ, :S2), ([-w / 2, -a / 2, 0.0], γ, :S3), ([w / 2, -a / 2, 0.0], γ, :S4)]
    east = [([a / 2, 0.0, t / 2], γ, :E1), ([a / 2, 0.0, -t / 2], γ, :E2), ([a / 2, w / 2, 0.0], γ, :E3), ([a / 2, -w / 2, 0.0], γ, :E4)]
    west = [([-a / 2, 0.0, t / 2], γ, :W1), ([-a / 2, 0.0, -t / 2], γ, :W2), ([-a / 2, w / 2, 0.0], γ, :W3), ([-a / 2, -w / 2, 0.0], γ, :W4)]
    excluder = [([0.0, 0.0, 0.0], γ, :central)]
    body_particles = Vector{Particle}()
    for (r, γγ, id) in [north; south; east; west; excluder]
        push!(body_particles, Particle(r, γγ, id))
    end
    temp_rigid_body = RigidBody(body_particles)
    rigid_body_pair = Vector{RigidBody}()
    for Δr = [box ./ 2 .- [a / 2, 0.0, 0.0], box ./ 2 .+ [a / 2, 0.0, 0.0]]
        new_rigid_body = deepcopy(temp_rigid_body)
        translate!(new_rigid_body, Δr)
        push!(rigid_body_pair, new_rigid_body)
    end
    all_particles = get_particle_list(rigid_body_pair)

    attractor_imatrix = create_interaction_matrix([[(Symbol("N$i"), Symbol("S$i")) for i = 1 : 4]; [(Symbol("E$i"), Symbol("W$i")) for i = 1 : 4]])
    attractors = get_particles_with_ids(all_particles, [Symbol("$(d)$(i)") for d in ["N", "S", "E", "W"], i = 1 : 4][:])
    attractor_bond_list = bind_closest(attractors, 0.01, attractor_imatrix)
    hb_rigid = HarmonicBond(k_corners, 0.0, attractor_bond_list, box; multithreaded = false)
    
    # Initialize integrator and system
    bodies = [lj_particle_pair ; m_particle_pair ; hb_particle_pair ; rigid_body_pair]
    brownian = Brownian(bodies, dt, kT, box; multithreaded = false)
    system = System(bodies, [lj, m, hb, hb_rigid], [lj_cell_list, m_cell_list], brownian)

    trajectories = Trajectories(save_interval)
    run_simulation!(system, trajectories, num_steps)
    #save!(system, "TEST_OUTPUT/system.out)
    save!(trajectories, "TEST_OUTPUT/trajectories_interaction_test.txt")
end
run!()