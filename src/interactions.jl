abstract type Interaction end

function wrap_displacement(Δx::Float64, Δy::Float64, Δz::Float64, box::SVector{3, Float64})

end

mutable struct LennardJones <: Interaction
    ϵ::Float64
    σ::Dict{Tuple{String, String}, Float64}
    r_cut::Dict{Tuple{String, String}, Float64}

    particles::Vector{Particle}
    neighbor_list::LinkedCellList
    box::SVector{3, Float64}
end

mutable struct Morse <: Interaction
    ϵ::Float64
    α::Dict{Tuple{String, String}, Float64}
    r_cut::Dict{Tuple{String, String}, Float64}

    particles::Vector{Particle}
    neighbor_list::LinkedCellList
    box::SVector{3, Float64}
end

mutable struct HarmonicBond <: Interaction
    k::Float64
    r0::SVector{3, Float64}

    bond_list::BondList
    box::SVector{3, Float64}
end