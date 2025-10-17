abstract type ExternalForce end

###########################################################################################################################################
# 
###########################################################################################################################################

struct ConstantForce <: ExternalForce
    force::SVector{3, Float64}

    particles::Vector{Particle}

    multithreaded::Bool
end
ConstantForce(force::Vector{Float64}, particles::Vector{Particle}; multithreaded::Bool = false) = ConstantForce(force, particles, multithreaded)

function compute_forces!(cf::ConstantForce)
    @use_threads cf.multithreaded for particle in cf.particles
        particle.force .+= cf.force
    end
end

###########################################################################################################################################
# 
###########################################################################################################################################

struct HarmonicTrap <: ExternalForce
    k::Float64
    r0::SVector{3, Int64}

    particles::Vector{Particle}
    
    multithreaded::Bool
end
HarmonicTrap(k::Float64, r0::Vector{Float64}, particles::Vector{Particle}; multithreaded::Bool = false) = HarmonicTrap(k, r0, particles, multithreaded)

function compute_forces!(ht::HarmonicTrap)
    @use_threads ht.multithreaded for particle in ht.particles
        particle.force[1] -= ht.k * (particle.position[1] - ht.r0[1])
        particle.force[2] -= ht.k * (particle.position[2] - ht.r0[2])
        particle.force[3] -= ht.k * (particle.position[3] - ht.r0[3])
    end
end