abstract type Integrator end

struct Brownian
    bodies::Vector{<:Body}

    dt::Float64
    kT::Float64

    box::SVector{3, Float64}

    multithreaded::Bool
end

function update_bodies!(integrator::Brownian)
    @use_threads integrator.multithreaded for body in integrator.bodies
        update_body!(body, integrator.dt, integrator.kT, box)
    end
end

function update_body!(body::RigidBody, dt::Float64, kT::Float64, box::SVector{3, Float64})
    fx, fy, fz = 0.0, 0.0, 0.0
    τx, τy, τz = 0.0, 0.0, 0.0
end

function update_body!(body::Particle, dt::Float64, kT::Float64, box::SVector{3, Float64})
    amp = sqrt(2 * kT / (particle.γ * dt))
    body.position[1] += body.force[1] / body.γ + amp * randn()
    body.position[2] += body.force[2] / body.γ + amp * randn()
    body.position[3] += body.force[3] / body.γ + amp * randn()

    if !(0.0 <= body.position[1] < box[1])
        body.position[1] = mod(body.position[1], box[1])
    end
    if !(0.0 <= body.position[2] < box[2])
        body.position[2] = mod(body.position[2], box[2])
    end
    if !(0.0 <= body.position[3] < box[3])
        body.position[3] = mod(body.position[3], box[3])
    end
end