abstract type NeighborList end

###########################################################################################################################################
# 
###########################################################################################################################################
struct BondList <: NeighborList
    bonds::Vector{Tuple{Particle, Vector{Particle}}}
end

function bind_closest(particles::Vector{Particle}, r_cut::Float64, interaction_matrix::DefaultDict{Tuple{Symbol, Symbol}, Bool, Bool})
    bonds = Vector{Tuple{Particle, Vector{Particle}}}()
    N_particles = length(particles)
    for i = 1 : N_particles
        particle_i = particles[i]
        temp_neighbors = Vector{Particle}()
        for j = 1 : N_particles
            particle_j = particles[j]
            if i != j && interaction_matrix[particle_i.id, particle_j.id] && sum((particle_i.position .- particle_j.position).^2) <= r_cut^2
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

mutable struct LinkedCellList <: NeighborList
    particles::Vector{Particle}
    start_index::Array{Int64, 3}
    next_index::Vector{Int64}

    cell_counts::SVector{3, Int64}
    cell_sizes::SVector{3, Float64}
    box::SVector{3, Float64}

    update_interval::Int64
    update_counter::Int64
end

function LinkedCellList(particles::Vector{Particle}, approx_cell_size::Float64, box::Vector{Float64}; 
                        update_interval::Int64 = 1, approx_padding::Float64 = 0.0)
    cell_counts = floor.(Int64, box ./ (approx_cell_size + approx_padding))
    cell_sizes = box ./ cell_counts

    start_index = -ones(Int64, cell_counts[1], cell_counts[2], cell_counts[3])
    next_index = -ones(Int64, length(particles))
    
    for (n, particle) in enumerate(particles)
        i, j, k = floor.(Int64, mod.(particle.position, box) ./ cell_sizes) .+ 1

        if start_index[i, j, k] > 0
            next_index[n] = start_index[i, j, k]
        end
        start_index[i, j, k] = n
    end

    return LinkedCellList(particles, start_index, next_index, cell_counts, cell_sizes, box, update_interval, 0), minimum(cell_sizes .- approx_cell_size)
end

function update_neighbor_list!(cell_list::LinkedCellList)
    if (cell_list.update_counter += 1) == cell_list.update_interval
        fill!(cell_list.start_index, -1)

        for (n, particle) in enumerate(cell_list.particles)
            i = floor(Int64, mod(particle.position[1], cell_list.box[1]) / cell_list.cell_sizes[1]) + 1
            j = floor(Int64, mod(particle.position[2], cell_list.box[2]) / cell_list.cell_sizes[2]) + 1
            k = floor(Int64, mod(particle.position[3], cell_list.box[3]) / cell_list.cell_sizes[3]) + 1

            if cell_list.start_index[i, j, k] > 0
                cell_list.next_index[n] = cell_list.start_index[i, j, k]
            else
                cell_list.next_index[n] = -1
            end
            cell_list.start_index[i, j, k] = n
        end
        cell_list.update_counter = 0
    end
end