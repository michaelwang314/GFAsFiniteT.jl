abstract type Body end

mutable struct Particle <: Body
    position::MVector{3, Float64}
    force::MVector{3, Float64}

    γ::Float64

    id::Int64
    body_id::Union{Int64, Nothing}
end

function Particle(position::V, γ::Float64, id::Int64, body_id::Union{Int64, Nothing}) where V <: AbstractVector
    return Particle(position, [0.0, 0.0, 0.0], γ, id, body_id)
end
Particle(position::V, id::Int64) where V <: AbstractVector = Particle(position, 1.0, id, nothing)

mutable struct RigidBody <: Body
    centroid::MVector{3, Float64}
    particles::Vector{Particle}
    
    axes::MVector{3, MVector{3, Float64}}
    γ_rot::SVector{3, Float64}
    γ_trans::Float64
end

function RigidBody(particles::Vector{Particle}, axes::Vector{Vector{Float64}}, γ_rot::Vector{Float64}, γ_trans::Float64)
    position = MVector{3}(0.0, 0.0, 0.0)
    for particle in particles
        position .+= particle.position
    end
    position ./= length(particles)

    return RigidBody(position, particles, axes, γ_rot, γ_trans)
end

function RigidBody(particles::Vector{Particle})
    centroid = MVector{3}(0.0, 0.0, 0.0)
    γ_total = 0.0
    I = [0 0 0; 0 0 0; 0 0 0]
    for particle in particles
        position = particle.position
        γ = particle.γ
        
        centroid .+= position * γ
        γ_total += γ
        I .+= γ * ((dot(position, position) * [1 0 0; 0 1 0; 0 0 1]) .- (position * transpose(position)))
    end
    centroid ./= γ_total

    γ_rot, evecs = eigvals(I), eigvecs(I)
    axes = [evecs[:, 1], evecs[:, 2], evecs[:, 3]]

    return RigidBody(centroid, particles, axes, γ_rot, γ_total)
end

function rotate!(body::RigidBody, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    sin, cos = sincos(θ)
    Rxx, Rxy, Rxz = cos + axis_x^2 * (1 - cos), axis_x * axis_y * (1 - cos) - axis_z * sin, axis_x * axis_z * (1 - cos) + axis_y * sin
    Ryx, Ryy, Ryz = axis_y * axis_x * (1 - cos) + axis_z * sin, cos + axis_y^2 * (1 - cos), axis_y * axis_z * (1 - cos) - axis_x * sin
    Rzx, Rzy, Rzz = axis_z * axis_x * (1 - cos) - axis_y * sin, axis_z * axis_y * (1 - cos) + axis_x * sin, cos + axis_z^2 * (1 - cos)

    for particle in body.particles
        rx, ry, rz = particle.position[1] - body.centroid[1], particle.position[2] - body.centroid[2], particle.position[3] - body.centroid[3]

        particle.position[1] = body.centroid[1] + Rxx * rx + Rxy * ry + Rxz * rz
        particle.position[2] = body.centroid[2] + Ryx * rx + Ryy * ry + Ryz * rz
        particle.position[3] = body.centroid[3] + Rzx * rx + Rzy * ry + Rzz * rz
    end

    for axis in body.axes
        axis[1], axis[2], axis[3] = Rxx * axis[1] + Rxy * axis[2] + Rxz * axis[3], Ryx * axis[1] + Ryy * axis[2] + Ryz * axis[3], Rzx * axis[1] + Rzy * axis[2] + Rzz * axis[3]
    end
end

function rotate!(body::RigidBody, origin_x::Float64, origin_y::Float64, origin_z::Float64, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    
end

function rotate!(bodies::Vector{RigidBody}, origin_x::Float64, origin_y::Float64, origin_z::Float64, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    for body in bodies
        rotate!(body, origin_x, origin_y, origin_z, axis_x, axis_y, axis_z, θ)
    end
end

function translate!(body::RigidBody, Δx::Float64, Δy::Float64, Δz::Float64)
    body.centroid[1] += Δx
    body.centroid[2] += Δy
    body.centroid[3] += Δz

    for particle in body.particles
        particle.position[1] += Δx
        particle.position[2] += Δy
        particle.position[3] += Δz
    end
end

function translate!(bodies::Vector{RigidBody}, Δx::Float64, Δy::Float64, Δz::Float64)
    for body in bodies
        translate!(body, Δx, Δy, Δz)
    end
end