abstract type ExternalForce end

struct ConstantForce
    force::SVector{3, Float64}

    particles::Vector{Particle}

    multithreaded::Bool
end
ConstantForce(force::Vector{Float64}, particles::Vector{Particle}) = ConstantForce(force, particles, false)

function compute_forces!(cf::ConstantForce)
    @use_threads cf.multithreaded for particle in cf.particles
        particle.force .+= cd.force
    end
end