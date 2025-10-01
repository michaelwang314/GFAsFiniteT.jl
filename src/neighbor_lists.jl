abstract type NeighborList end

###########################################################################################################################################
# 
###########################################################################################################################################
struct BondList <: NeighborList
    bonds::Vector{Tuple{Particle, Vector{Particle}}}
end

function bind_closest(particles::Vector{Particle}, r_cut::Float64, interaction_matrix::Dict{Tuple{String, String}, Bool})
    bonds = Vector{Tuple{Particle, Vector{Particle}}}()
    N_particles = length(particles)
    for i = 1 : N_particles
        particle_i = particles[i]
        temp_neighbors = Vector{Particle}()
        for j = 1 : N_particles
            particle_j = particles[j]
            if i!=j && get(interaction_matrix, (particle_i.id, particle_j.id), false) && sum((particle_i.position .- particle_j.position).^2) <= r_cut^2
                push!(temp_neighbors, particle_j)
            end
        end

        if !isempty(temp_neighbors)
            push!(bonds, (particle_i, temp_neighbors))
        end
    end

    return BondList(bonds)
end

###########################################################################################################################################
# 
###########################################################################################################################################

struct LinkedCellList <: NeighborList
    particles::Vector{Particle}
    start_index::Array{Int64, 3}
    next_index::Vector{Int64}

    cell_counts::SVector{3, Int64}
    cell_spacings::SVector{3, Float64}
    box::SVector{3, Float64}
end

function LinkedCellList(particles::Vector{Particle}, approx_cell_spacings::Vector{Float64}, box::Vector{Float64})
    cell_counts = floor.(Int64, box ./ approx_cell_spacings)
    cell_spacings = box ./ cell_counts

    start_index = -ones(Int64, cell_counts[1], cell_counts[2], cell_counts[3])
    next_index = -ones(Int64, length(particles))
    
    for (n, particle) in enumerate(particles)
        i, j, k = floor.(Int64, mod.(particle.position, box) ./ cell_spacings)

        if start_index[i, j, k] > 0
            next_index[n] = start_index[i, j, k]
        end
        start_index[i, j, k] = n
    end

    return LinkedCellList(particles, start_index, next_index, cell_counts, cell_spacings, box)
end

function update_cell_list!(cell_list::LinkedCellList)
    fill!(cell_list.start_index, -1)
    fill!(cell_list.next_index, -1)

    for (n, particle) in enumerate(cell_list.particles)
        i = floor(Int64, mod(particle.position[1], cell_list.box[1]) / cell_list.cell_spacings[1]) + 1
        j = floor(Int64, mod(particle.position[2], cell_list.box[2]) / cell_list.cell_spacings[2]) + 1
        k = floor(Int64, mod(particle.position[3], cell_list.box[3]) / cell_list.cell_spacings[3]) + 1

        if cell_list.start_index[i, j, k] > 0
            cell_list.next_index[n] = cell_list.start_index[i, j, k]
        end
        cell_list.start_index[i, j, k] = n
    end
end