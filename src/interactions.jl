abstract type Interaction end

function wrap_displacement(Δx::Float64, Δy::Float64, Δz::Float64, box::SVector{3, Float64})

end

mutable struct LennardJones <: Interaction
    ϵ::Float64
    σ::Float64
    r_cut::Float64

    particles::Vector{Particle}
    neighbor_list::LinkedCellList
end

mutable struct Morse <: Interaction
    ϵ::Float64
    α::Float64
    r_cut::Float64

    particles::Vector{Particle}
    neighbor_list::LinkedCellList
end

mutable struct HarmonicBond <: Interaction
    k::Float64
    r0::SVector{3, Float64}

    bond_list::BondList
end