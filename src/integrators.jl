abstract type Integrator end

struct Brownian <: Integrator
    bodies::Vector{<:Body}

    dt::Float64
    kT::Float64

    box::SVector{3, Float64}

    multithreaded::Bool
end

function update_bodies!(integrator::Brownian)
    @use_threads integrator.multithreaded for body in integrator.bodies
        update_body!(body, integrator.dt, integrator.kT, integrator.box)
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
    τ1_scaled = (τx_total * body.axes[1][1] + τy_total * body.axes[1][2] + τz_total * body.axes[1][3]) / body.moments[1]
    τ2_scaled = (τx_total * body.axes[2][1] + τy_total * body.axes[2][2] + τz_total * body.axes[2][3]) / body.moments[2]
    τ3_scaled = (τx_total * body.axes[3][1] + τy_total * body.axes[3][2] + τz_total * body.axes[3][3]) / body.moments[3]
    ωx = τ1_scaled * body.axes[1][1] + τ2_scaled * body.axes[2][1] + τ3_scaled * body.axes[3][1]
    ωy = τ1_scaled * body.axes[1][2] + τ2_scaled * body.axes[2][2] + τ3_scaled * body.axes[3][2]
    ωz = τ1_scaled * body.axes[1][3] + τ2_scaled * body.axes[2][3] + τ3_scaled * body.axes[3][3]
    ω = sqrt(ωx^2 + ωy^2 + ωz^2)

    dt_scaled = dt / body.γ_trans
    translate!(body, dt_scaled * fx_total, dt_scaled * fy_total, dt_scaled * fz_total)
    if ω > 0.0
        rotate!(body, ωx / ω, ωy / ω, ωz / ω, dt * ω)
    end

    Δx, Δy, Δz = 0.0, 0.0, 0.0
    if !(0.0 <= body.centroid[1] < box[1])
        Δx = mod(body.centroid[1], box[1]) - body.centroid[1]
    end
    if !(0.0 <= body.centroid[2] < box[2])
        Δy = mod(body.centroid[2], box[2]) - body.centroid[2]
    end
    if !(0.0 <= body.centroid[3] < box[3])
        Δz = mod(body.centroid[3], box[3]) - body.centroid[3]
    end
    if Δx != 0.0 || Δy != 0.0 || Δz != 0.0
        translate!(body, Δx, Δy, Δz)
    end
end

function update_body!(body::Particle, dt::Float64, kT::Float64, box::SVector{3, Float64})
    amp = sqrt(2 * kT / (body.γ * dt))
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