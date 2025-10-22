struct System
    bodies::Vector{<:Body}
    interactions::Vector{<:Interaction}
    neighbor_lists::Vector{<:NeighborList}
    external_forces::Vector{<:ExternalForce}
    integrator::Integrator
end
System(bodies::Vector{<:Body}, interactions::Vector{<:Interaction}, neighbor_lists::Vector{<:NeighborList}, integrator::Integrator) = System(bodies, interactions, neighbor_lists, ExternalForce[], integrator)

mutable struct Trajectories
    history::Vector{Vector{<:Body}}

    start::Int64
    period::Int64
end
Trajectories(period::Int64; start::Int64 = 1) = Trajectories(Vector{Vector{<:Body}}(), start, period)

function run_simulation!(system::System, trajectories::Union{Trajectories, Nothing}, num_steps::Int64; message_interval::Union{Float64, Nothing} = 10.0)
    prev_step = 0
    time_elapsed = 0.0
    interval_start = time()
    for step = 1 : num_steps
        for interaction in system.interactions
            compute_forces!(interaction)
        end

        for external_force in system.external_forces
            compute_forces!(external_force)
        end

        update_bodies!(system.integrator)
        
        for neighbor_list in system.neighbor_lists
            update_neighbor_list!(neighbor_list)
        end

        if !isnothing(trajectories) && (step - trajectories.start) % trajectories.period == 0
            push!(trajectories.history, deepcopy(system.bodies))
        end

        interval_time = time() - interval_start
        if !isnothing(message_interval) && (interval_time > message_interval || step == num_steps)
            time_elapsed += interval_time
            rate = (step - prev_step) / interval_time
            println(hr_min_sec(time_elapsed), " | ",
                    step, "/", num_steps, " (", round(step / num_steps * 100, digits = 1), "%) | ",
                    round(rate, digits = 1), " steps/s | ",
                    hr_min_sec((num_steps - step) / rate))
            prev_step = step
            interval_start = time()
        end
    end
    println("Average steps/s: ", round(num_steps / time_elapsed, digits = 1))
end
run_simulation!(system::System, num_steps::Int64; message_interval::Union{Float64, Nothing} = 10.0) = run_simulation!(system, nothing, num_steps; message_interval = message_interval)

function initialize_random_positions(N::Int64, r_cut::Float64, box::Vector{Float64}; max_attempts::Int64 = 100, excluded_region::Union{Function, Nothing} = nothing)
    positions = Vector{Vector{Float64}}()

    cell_counts = floor.(Int64, box ./ r_cut)
    cell_sizes = box ./ cell_counts
    cell_ids = DefaultDict{Tuple{Int64, Int64, Int64}, Vector{Int64}}(Vector{Int64}())

    n = 0
    attempts = 0
    while n < N && attempts < max_attempts 
        x_new, y_new, z_new = rand(3) .* box

        if !isnothing(excluded_region) && excluded_region(x_new, y_new, z_new)
           attempts += 1
           continue
        end

        is_overlapping = false
        i, j, k = floor.(Int64, [x_new, y_new, z_new] ./ cell_sizes) .+ 1
        for Δi = -1 : 1, Δj = -1 : 1, Δk = -1 : 1
            iΔi, jΔj, kΔk = mod.([i + Δi, j + Δj, k + Δk], cell_counts) .+ 1
            for check_id in cell_ids[iΔi, jΔj, kΔk]
                x_check, y_check, z_check = positions[check_id]

                if (x_check - x_new)^2 + (y_check - y_new)^2 + (z_check - z_new)^2 < r_cut^2
                    is_overlapping = true
                    break
                end
            end
            if is_overlapping
                break
            end
        end
        if is_overlapping
            attempts += 1
            continue
        end

        n += 1
        push!(cell_ids[i, j, k], n)
        push!(positions, [x_new, y_new, z_new])
    end

    if attempts == max_attempts
        println("Maximum attempts reached when initializing random positions.  Not all positions were generated.")
    end

    return positions
end

function initialize_lattice_positions()
    # to be finished
end

###########################################################################################################################################
# 
###########################################################################################################################################

function hr_min_sec(time::Float64)
    hours = trunc(Int64, time / 3600.0)
    minutes = trunc(Int64, mod(time, 3600.0) / 60.0)
    seconds = trunc(Int64, mod(time, 60.0))

    return string(hours < 10 ? "0" : "", hours, 
                  minutes < 10 ? ":0" : ":", minutes, 
                  seconds < 10 ? ":0" : ":", seconds)
end

function get_ext(filename::String)
    return filename[findlast(isequal('.'), filename) : end]
end

function save!(system::System, filename::String)
    if !isdir(dirname(filename))
        mkpath(dirname(filename))
    end

    open(filename, "w") do io
        serialize(io, system)
    end

    println("$(filename) saved!")
end

function save!(trajectories::Trajectories, filename::String)
    if !isdir(dirname(filename))
        mkpath(dirname(filename))
    end

    if get_ext(filename) == ".txt"
        open(filename, "w") do io
            write(io, "t, x, y, z, id, body id\n")
            for (t, bodies) in enumerate(trajectories.history)
                for particle in get_particle_list(bodies)
                    write(io, "$(t), $(particle.position[1]), $(particle.position[2]), $(particle.position[3]), $(particle.id), $(particle.body_id)\n")
                end
            end
        end
    else
        open(filename, "w") do io
            serialize(io, trajectories)
        end
    end
    println("$(filename) saved!")
end

function load(filename::String)
    object = begin
        open(filename, "r") do io
            deserialize(io)
        end
    end
    println("$(filename) loaded!")

    return object
end
