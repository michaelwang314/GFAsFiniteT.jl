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