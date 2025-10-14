using GFAsFiniteT

function run!()
    box = [20.0, 20.0, 20.0]
    num_steps = 1000000
    save_interval = trunc(Int64, num_steps / 10000)
    dt = 0.0001
    kT = 0.1

    ϵ = 1.0
    σ = 1.0
    lj_r_cut = 2 * 2^(1 / 6) * σ
    lj_particle_pair = [Particle(box ./ 2 .- [lj_r_cut / 4, 0.0, 0.0], "lj particle"), Particle(box ./ 2 .+ [lj_r_cut / 4, 0.0, 0.0], "lj particle")]
    lj_cell_list = LinkedCellList(lj_particle_pair, lj_r_cut, box)
    lj_params = [("lj particle", "lj particle", Dict("ϵ" => ϵ, "σ" => σ, "r_cut" => lj_r_cut))]
    lj = LennardJones(lj_params, lj_particle_pair, lj_cell_list, box, false)

    D0 = 5.0
    α = 10.0
    r0 = 0.1
    m_r_cut = 2.0
    m_particle_pair = [Particle(box ./ 2, "morse particle"), Particle(box ./ 2, "morse particle")]
    m_cell_list = LinkedCellList(m_particle_pair, m_r_cut, box)
    m_params = [("morse particle", "morse particle", Dict("D0" => D0, "α" => α, "r0" => r0, "r_cut" => m_r_cut))]
    m = Morse(m_params, m_particle_pair, m_cell_list, box, false)

    k = 10.0
    r0 = 0.1
    hb_particle_pair = [Particle(box ./ 2 .- [r0 / 2, 0.0, 0.0], "harmonic particle"), Particle(box ./ 2 .+ [r0 / 2, 0.0, 0.0], "harmonic particle")]
    hb_imatrix = create_interaction_matrix([("harmonic particle", "harmonic particle")])
    bond_list = bind_closest(hb_particle_pair, 2 * r0, hb_imatrix)
    hb = HarmonicBond(k, r0, bond_list, box, false)

    bodies = [lj_particle_pair ; m_particle_pair ; hb_particle_pair]
    brownian = Brownian(bodies, dt, kT, box, false)
    system = System(bodies, [lj, m, hb], [lj_cell_list, m_cell_list], brownian)

    trajectories = Trajectories(save_interval)
    run_simulation!(system, trajectories, num_steps)
    #save_system!(system, "TEST_OUTPUT/system.out")
    export_trajectories!(trajectories, "TEST_OUTPUT/trajectories_interaction_test.txt")
end
run!()