struct System
    bodies::Vector{<:Body}
    interactions::Vector{<:Interaction}
    external_forces::Vector{<:ExternalForce}
    integrator::Integrator

    box::SVector{3, Float64}
end

mutable struct Trajectories
    history::Vector{Vector{<:Body}}

    start::Int64
    period::Int64
end
Trajectories(start::Int64, period::Int64) = Trajectories(Vector{Vector{<:Body}}(), start, period)

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
    time_elasped = 0.0
    interval_start = time()
    for step = 1 : num_steps
        for interaction in system.interactions
            compute_forces!(interactions)
        end
        for external_force in system.external_forces
            compute_forces!(external_force)
        end
        update_bodies!(system.bodies)

        if !isnothing(trajectories) && (step - trajectories.start) % trajectories.period
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