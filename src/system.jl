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
Trajectories(start::Int64, period::Int64) = Trajectories(Vector{Vector{<:Body}}(), start, period)
Trajectories(period::Int64) = Trajectories(1, period)

function hr_min_sec(time::Float64)
    hours = trunc(Int64, time / 3600.0)
    minutes = trunc(Int64, mod(time, 3600.0) / 60.0)
    seconds = trunc(Int64, mod(time, 60.0))

    return string(hours < 10 ? "0" : "", hours, 
                  minutes < 10 ? ":0" : ":", minutes, 
                  seconds < 10 ? ":0" : ":", seconds)
end

function run_simulation!(system::System, trajectories::Union{Trajectories, Nothing}, num_steps::Int64; message_interval::Float64 = 10.0)
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
        if interval_time > message_interval || step == num_steps
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
end
run_simulation!(system::System, num_steps::Int64; message_interval::Float64 = 10.0) = run_simulation!(system, nothing, num_steps; message_interval = message_interval)

function save_system!(system::System, filename::String)
    if !isdir(dirname(filename))
        mkpath(dirname(filename))
    end

    open(filename, "w") do io
        serialize(io, system)
    end
end

function load_system(filename::String)
    return begin
        open(filename, "r") do io
            deserialize(io)
        end
    end
end

function export_trajectories!(trajectories::Trajectories, filename::String, multi_particle_types::Vector{DataType})
    if !isdir(dirname(filename))
        mkpath(dirname(filename))
    end

    open(filename, "w") do io
        write(io, "t, x, y, z, id, body id\n")
        for (t, bodies) in enumerate(trajectories.history)
            all_particles = get_particle_list(bodies, multi_particle_types)

            for particle in all_particles
                write(io, "$(t), $(particle.position[1]), $(particle.position[2]), $(particle.position[3]), $(particle.id), $(particle.body_id)\n")
            end
        end
    end
end
