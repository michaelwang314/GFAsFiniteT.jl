abstract type NeighborList end

struct BondList <: NeighborList
    bonds::Vector{Tuple{Particle, Vector{Particle}}}
end

struct LinkedCellList <: NeighborList
    particles::Vector{Particle}
    start_id::Array{Int64, 3}
    next_id::Vector{Int64}

    cell_count::SVector{Int64, 3}
    cell_spacing::SVector{Float64, 3}
end

function create_interaction_matrix(pairs::Vector{Tuple{String, String}})
    interaction_matrix = Dict{Tuple{String, String}, Bool}()
    for (id1, id2) in pairs
        interaction_matrix[(id1, id2)] = interaction_matrix[(id2, id1)] = true
    end

    return interaction_matrix
end

function bind_closest(particles::Vector{Particle}, r_cut::Float64, interaction_matrix::Dict{Tuple{String, String}, Bool})
    bonds = Vector{Tuple{Particle, Vector{Particle}}}()
    N_particles = length(particles)
    for i = 1 : N_particles
        particle_i = particles[i]
        temp_neighbors = Vector{Particles}()
        for j = 1 : N_particles
            particle_j = particles[j]
            if i != j && interaction_matrix[(particle_i.id, particle_j.id)] && sum((particle_i.position .- particle_j.position).^2) <= r_cut^2
                push!(temp_neighbors, particle_j)
            end
        end

        if !isempty(temp_neighbors)
            push!(bonds, temp_neighbors)
        end
    end

    return BondList(bonds)
end

function rigid_bodies_to_particle_list(bodies::Vector{RigidBody})
    particles = Vector{Particles}()
    for body in bodies
        append!(particles, body.particles)
    end
    return particles
end

function get_particles_with_id(particles::Vector{Particle}, ids::Vector{String})
    subset = Vector{Particle}()
    for particle in particles
        if particle.id in ids
            push!(subset, particle)
        end
    end
    return subset
end

