using GFAsFiniteT
#=using ArgParse

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--w"
            arg_type = Float64
            default = 1.0
        "--numlinkers"
            arg_type = Int64
            default = 0
        "--traj_filename"
            arg_type = String
            default = "TEST_OUTPUTS/trajectories_tetramer_test.txt"
    end

    return parse_args(s)
end=#

function run!()
    #=cmd_args = parse_cmd()
    w_scale = cmd_args["w"]
    num_linkers = cmd_args["numlinkers"]
    filename = cmd_args["traj_filename"]=#
    w_scale = 0.0
    num_linkers = 0
    sys_filename = "TEST_OUTPUTS/system_tetramer_test.out"
    traj_filename = "TEST_OUTPUTS/trajectories_tetramer_test.txt"
    
    box = [5.0, 5.0, 5.0]
    num_steps = 1000000
    save_interval = trunc(Int64, num_steps / 10000)
    dt = 0.0005
    kT = 0.1

    γ_linker = 1.0;
    γ_constituents = 1.0;

    a = 1.0
    t = a / sqrt(2)
    w = w_scale * t
    k_corners = 20.0

    r_linker_site = 0.125
    r_linker = 0.125
    r_excluder = a / (2 * sqrt(2))

    north = [([0.0, a / 2, t / 2], γ_constituents, :N1), ([0.0, a / 2, -t / 2], γ_constituents, :N2), ([-w / 2, a / 2, 0.0], γ_constituents, :N3), ([w / 2, a / 2, 0.0], γ_constituents, :N4)]
    south = [([0.0, -a / 2, t / 2], γ_constituents, :S1), ([0.0, -a / 2, -t / 2], γ_constituents, :S2), ([-w / 2, -a / 2, 0.0], γ_constituents, :S3), ([w / 2, -a / 2, 0.0], γ_constituents, :S4)]
    east = [([a / 2, 0.0, t / 2], γ_constituents, :E1), ([a / 2, 0.0, -t / 2], γ_constituents, :E2), ([a / 2, w / 2, 0.0], γ_constituents, :E3), ([a / 2, -w / 2, 0.0], γ_constituents, :E4)]
    west = [([-a / 2, 0.0, t / 2], γ_constituents, :W1), ([-a / 2, 0.0, -t / 2], γ_constituents, :W2), ([-a / 2, w / 2, 0.0], γ_constituents, :W3), ([-a / 2, -w / 2, 0.0], γ_constituents, :W4)]
    linker_sites = [([a / 4, a / 4, 0.0], γ_constituents, :linker_site), ([-a / 4, a / 4, 0.0], γ_constituents, :linker_site), ([-a / 4, -a / 4, 0.0], γ_constituents, :linker_site), ([a / 4, -a / 4, 0.0], γ_constituents, :linker_site)]
    excluders = [([0.0, 0.0, 0.0], γ_constituents, :central)]
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
        set_body_ids!(new_rigid_body, Symbol("body_$(i += 1)"))
        push!(bodies, new_rigid_body)
    end

    excluded(x, y, z) = box[1] / 2 - 1.0 < x < box[1] / 2 + 1.0 && box[2] / 2 - 1.0 < y < box[2] / 2 + 1.0 && box[3] / 2 - 0.5 < z < box[3] / 2 + 0.5
    for position in initialize_random_positions(num_linkers, 2^(1 / 6) * 2 * r_linker, box; excluded_region = excluded)
        push!(bodies, Particle(position, γ_linker, :linker))
    end
    
    all_particles = get_particle_list(bodies)
    
    #=non_hb_particles = get_particles_with_ids(all_particles, [:central, :linker_site, :linker])
    max_r_cut = maximum([2^(1 / 6) * 2 * r_excluder, 2^(1 / 6) * (r_excluder + r_linker), 2^(1 / 6) * 2 * r_linker, r_linker_site + r_linker])
    cell_list, padding = LinkedCellList(non_hb_particles, max_r_cut, box)
    cell_list.update_interval = maximum([1, floor(Int64, 0.75 * padding^2 / (2 * 0.1 * 0.0001 / 1.0))])
    println("CellList update interval set to $(cell_list.update_interval)")

    m_params = [(:linker_site, :linker, Dict(:D0 => 0.5, :α => 2 / (r_linker_site + r_linker), :r0 => 0.0, :r_cut => r_linker_site + r_linker))]
    m = Morse(m_params, cell_list, box; multithreaded = false)

    #(:central, :linker, Dict(:ϵ => 1.0, :σ => r_excluder + r_linker, :r_cut => 2^(1 / 6) * (r_linker + r_excluder))),
    lj_params = [(:central, :central, Dict(:ϵ => 1.0, :σ => 2 * r_excluder, :r_cut => 2^(1 / 6) * 2 * r_excluder)),
                 (:linker, :linker, Dict(:ϵ => 1.0, :σ => 2 * r_linker, :r_cut => 2^(1 / 6) * 2 * r_linker))]
    lj = LennardJones(lj_params, cell_list, box; multithreaded = true)=#
    
    small_particles = get_particles_with_ids(all_particles, [:linker_site, :linker])
    small_cell_list, padding = LinkedCellList(small_particles, maximum([r_linker_site + r_linker, 2^(1 / 6) * 2 * r_linker]), box)
    small_cell_list.update_interval = maximum([1, floor(Int64, 0.75 * padding^2 / (2 * kT * dt / 1.0))])
    println("CellList update interval set to $(small_cell_list.update_interval)")

    large_particles = get_particles_with_ids(all_particles, [:central])
    large_cell_list, padding = LinkedCellList(large_particles, 2^(1 / 6) * 2 * r_excluder, box)
    large_cell_list.update_interval = maximum([1, floor(Int64, 0.75 * padding^2 / (2 * kT * dt / 10.0))])
    println("CellList update interval set to $(large_cell_list.update_interval)")

    m_params = [(:linker_site, :linker, Dict(:D0 => 0.5, :α => 2 / (r_linker_site + r_linker), :r0 => 0.0, :r_cut => r_linker_site + r_linker))]
    m = Morse(m_params, small_cell_list, box; multithreaded = false)

    linker_lj_params = [(:linker, :linker, Dict(:ϵ => 1.0, :σ => 2 * r_linker, :r_cut => 2^(1 / 6) * 2 * r_linker))]
    linker_lj = LennardJones(linker_lj_params, small_cell_list, box; multithreaded = true)

    excluder_lj_params = [(:central, :central, Dict(:ϵ => 1.0, :σ => 2 * r_excluder, :r_cut => 2^(1 / 6) * 2 * r_excluder))]
    excluder_lj = LennardJones(excluder_lj_params, large_cell_list, box; multithreaded = false)

    attractors = get_particles_with_ids(all_particles, [Symbol("$(d)$(i)") for d in ["N", "S", "E", "W"], i = 1 : 4][:])
    attractor_imatrix = create_interaction_matrix([[(Symbol("N$i"), Symbol("S$i")) for i = 1 : 4]; [(Symbol("E$i"), Symbol("W$i")) for i = 1 : 4]])
    hb_bond_list = bind_closest(attractors, 0.01, attractor_imatrix)
    hb = HarmonicBond(k_corners, 0.0, hb_bond_list, box; multithreaded = false)

    brownian = Brownian(bodies, dt, kT, box; multithreaded = true)

    #system = System(bodies, [lj, hb, m], [cell_list], brownian)
    system = System(bodies, [excluder_lj, linker_lj, hb, m], [small_cell_list, large_cell_list], brownian)
    trajectories = Trajectories(save_interval)
    run_simulation!(system, trajectories, num_steps)
    save!(system, sys_filename)
    save!(trajectories, traj_filename)
end
run!()