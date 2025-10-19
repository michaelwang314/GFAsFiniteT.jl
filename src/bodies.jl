abstract type Body end

###########################################################################################################################################
# Particle
###########################################################################################################################################

mutable struct Particle <: Body
    position::MVector{3, Float64}
    force::MVector{3, Float64}

    γ::Float64

    id::Symbol
    body_id::Symbol
end
Particle(position::Vector{Float64}, γ::Float64, id::Symbol) = Particle(position, [0.0, 0.0, 0.0], γ, id, :none)

###########################################################################################################################################
# Rigid body
###########################################################################################################################################
"""
`RigidBody` is a rigid collection of `Particles`

Requires a list of particles.  Centroid and moment of friction/inertia computed from particle positions. 
"""
mutable struct RigidBody <: Body
    centroid::MVector{3, Float64}
    particles::Vector{Particle}
    
    axes::MVector{3, MVector{3, Float64}}
    moments::SVector{3, Float64}
    γ_trans::Float64

    id::Symbol

    function RigidBody(centroid::Vector{Float64}, particles::Vector{Particle}, axes::Vector{Vector{Float64}}, moments::Vector{Float64}, γ_trans::Float64, id::Symbol)
        for particle in particles
            particle.body_id = id
        end
    
        return new(centroid, particles, axes, moments, γ_trans, id)
    end
end

function RigidBody(particles::Vector{Particle}, id::Symbol)
    centroid = zeros(3)
    γ_total = 0.0
    I = zeros(3, 3)
    for particle in particles
        position = particle.position
        γ = particle.γ
        
        centroid .+= position * γ
        γ_total += γ
        I .+= γ * ((dot(position, position) * diagm([1.0, 1.0, 1.0])) .- (position * transpose(position)))
    end
    centroid ./= γ_total

    moments, evecs = eigvals(I), eigvecs(I)
    axes = [evecs[:, 1], evecs[:, 2], evecs[:, 3]]

    return RigidBody(centroid, particles, axes, moments, γ_total, id)
end
RigidBody(particles::Vector{Particle}) = RigidBody(particles, :none)

"""
    rotate!(body, axis_x, axis_y, axis_z, θ)

Rotate `body` about an axis `[axis_x, axis_y, axis_z]` around its centroid by an angle `θ`
"""
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
rotate!(body::RigidBody, axis::Vector{Float64}, θ::Float64) = rotate!(body, axis[1], axis[2], axis[3], θ)

"""
    rotate!(body, origin_x, origin_y, origin_z, axis_x, axis_y, axis_z, θ)

Rotate a `body` about an axis `[axis_x, axis_y, axis_z]` around an origin `[origin_x, origin_y, origin_z]` by an angle `θ` 
"""
function rotate!(body::RigidBody, origin_x::Float64, origin_y::Float64, origin_z::Float64, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    sin, cos = sincos(θ)
    Rxx, Rxy, Rxz = cos + axis_x^2 * (1 - cos), axis_x * axis_y * (1 - cos) - axis_z * sin, axis_x * axis_z * (1 - cos) + axis_y * sin
    Ryx, Ryy, Ryz = axis_y * axis_x * (1 - cos) + axis_z * sin, cos + axis_y^2 * (1 - cos), axis_y * axis_z * (1 - cos) - axis_x * sin
    Rzx, Rzy, Rzz = axis_z * axis_x * (1 - cos) - axis_y * sin, axis_z * axis_y * (1 - cos) + axis_x * sin, cos + axis_z^2 * (1 - cos)

    rx, ry, rz = body.centroid[1] - origin_x, body.centroid[2] - origin_y, body.centroid[3] - origin_z

    body.centroid[1] = origin_x + Rxx * rx + Rxy * ry + Rxz * rz
    body.centroid[2] = origin_y + Ryx * rx + Ryy * ry + Ryz * rz
    body.centroid[3] = origin_z + Rzx * rx + Rzy * ry + Rzz * rz

    for particle in body.particles
        rx, ry, rz = particle.position[1] - origin_x, particle.position[2] - origin_y, particle.position[3] - origin_z

        particle.position[1] = origin_x + Rxx * rx + Rxy * ry + Rxz * rz
        particle.position[2] = origin_y + Ryx * rx + Ryy * ry + Ryz * rz
        particle.position[3] = origin_z + Rzx * rx + Rzy * ry + Rzz * rz
    end

    for axis in body.axes
        axis[1], axis[2], axis[3] = Rxx * axis[1] + Rxy * axis[2] + Rxz * axis[3], Ryx * axis[1] + Ryy * axis[2] + Ryz * axis[3], Rzx * axis[1] + Rzy * axis[2] + Rzz * axis[3]
    end
end
rotate!(body::RigidBody, origin::Vector{Float64}, axis::Vector{Float64}, θ) = rotate!(body, origin[1], origin[2], origin[3], axis[1], axis[2], axis[3], θ)

"""
    rotate!(bodies, ...)

Rotate multiple bodies collectively
"""
function rotate!(bodies::Vector{RigidBody}, origin_x::Float64, origin_y::Float64, origin_z::Float64, axis_x::Float64, axis_y::Float64, axis_z::Float64, θ::Float64)
    for body in bodies
        rotate!(body, origin_x, origin_y, origin_z, axis_x, axis_y, axis_z, θ)
    end
end
rotate!(bodies::Vector{RigidBody}, origin::Vector{Float64}, axis::Vector{Float64}, θ) = rotate!(bodies, origin[1], origin[2], origin[3], axis[1], axis[2], axis[3], θ)

"""
    translate!(body, Δx, Δy, Δz)

Translate a `body` by `[Δx, Δy, Δz]`
"""
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
translate!(body::RigidBody, Δr::Vector{Float64}) = translate!(body, Δr[1], Δr[2], Δr[3])

"""
    translate!(bodies, ...)

Translate multiple bodies collectively
"""
function translate!(bodies::Vector{RigidBody}, Δx::Float64, Δy::Float64, Δz::Float64)
    for body in bodies
        translate!(body, Δx, Δy, Δz)
    end
end
translate!(bodies::Vector{RigidBody}, Δr::Vector{Float64}) = translate!(bodies, Δr[1], Δr[2], Δr[3])

###########################################################################################################################################
# Additional functions
###########################################################################################################################################
"""
    set_body_ids!(body, id)

Set `id` of `body` as well as `body_id` of constituent particles
"""
function set_body_ids!(body::RigidBody, id::Symbol)
    body.id = id
    for particle in body.particles
        particle.body_id = id
    end
end

function get_particle_list(bodies::Vector{<:Body})
    particles = Vector{Particle}()
    for body in bodies
        if hasfield(typeof(body), :particles)
            append!(particles, body.particles)
        else
            push!(particles, body)
        end
    end
    return particles
end

function get_particles_with_ids(particles::Vector{Particle}, ids::Vector{Symbol})
    subset = Vector{Particle}()
    for particle in particles
        if particle.id in ids
            push!(subset, particle)
        end
    end
    return subset
end