abstract type Interaction end

###########################################################################################################################################
# 
###########################################################################################################################################

mutable struct LennardJones <: Interaction
    ϵ::Float64
    σ::Dict{Tuple{String, String}, Float64}
    r_cut::Dict{Tuple{String, String}, Float64}

    neighbor_info::Tuple{Vector{Particle}, LinkedCellList, Dict{Tuple{String, String}, Bool}}
    box::SVector{3, Float64}
end

function compute_force!(lj::LennardJones)
    
end

###########################################################################################################################################
# 
###########################################################################################################################################

mutable struct Morse <: Interaction
    ϵ::Float64
    α::Dict{Tuple{String, String}, Float64}
    r_cut::Dict{Tuple{String, String}, Float64}

    neighbor_info::Tuple{Vector{Particle}, LinkedCellList, Dict{Tuple{String, String}, Bool}}
    box::SVector{3, Float64}
end

function compute_force!(m::Morse)

end

###########################################################################################################################################
# 
###########################################################################################################################################

mutable struct HarmonicBond <: Interaction
    k::Float64
    r0::SVector{3, Float64}

    bond_list::BondList
    box::SVector{3, Float64}
end

function compute_force!(hb::HarmonicBond)
    
end

###########################################################################################################################################
# Additional functions
###########################################################################################################################################

function wrap_displacement(Δx::Float64, Δy::Float64, Δz::Float64, box::SVector{3, Float64})

end

function create_interaction_matrix(pairs::Vector{Tuple{String, String}})
    interaction_matrix = Dict{Tuple{String, String}, Bool}()
    for (id1, id2) in pairs
        interaction_matrix[(id1, id2)] = interaction_matrix[(id2, id1)] = true
    end

    return interaction_matrix
end