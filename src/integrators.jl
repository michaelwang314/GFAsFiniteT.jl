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
    fx_total, fy_total, fz_total = 0.0, 0.0, 0.0
    τx_total, τy_total, τz_total = 0.0, 0.0, 0.0
    for particle in body.particles
        amp = sqrt(2 * kT * particle.γ / dt)
        fx_total += (fx = particle.force[1] + amp * randn())
        fy_total += (fy = particle.force[2] + amp * randn())
        fz_total += (fz = particle.force[3] + amp * randn())

        rx, ry, rz = particle.position[1] - body.centroid[1], particle.position[2] - body.centroid[2], particle.position[3] - body.centroid[3]
        τx_total += ry * fz - fy * rz 
        τy_total += -(rx * fz - fx * rz) 
        τz_total += rx * fy - fx * ry

        fill!(particle.force, 0.0)
    end

    translate!(body, dt * fx_total, dt * fy_total, dt * fz_total)
    τ1 = τx_total * body.axes[1][1] + τy_total * body.axes[1][2] + τz_total * body.axes[1][3]
    τ2 = τx_total * body.axes[2][1] + τy_total * body.axes[2][2] + τz_total * body.axes[2][3]
    τ3 = τx_total * body.axes[3][1] + τy_total * body.axes[3][2] + τz_total * body.axes[3][3]
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