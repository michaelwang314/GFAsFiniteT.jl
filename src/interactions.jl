abstract type Interaction end

###########################################################################################################################################
# 
###########################################################################################################################################

mutable struct LennardJones <: Interaction
    ϵ::Dict{Tuple{String, String}, Float64}
    σ::Dict{Tuple{String, String}, Float64}
    r_cut::Dict{Tuple{String, String}, Float64}

    particles::Vector{Particle}
    neighbor_list::LinkedCellList
    interaction_matrix::DefaultDict{Tuple{String, String}, Bool, Bool}
    box::SVector{3, Float64}

    multithreaded::Bool
end

function LennardJones(params::Vector{Tuple{String, String, Dict{String, Float64}}}, particles::Vector{Particle}, neighbor_list::LinkedCellList, 
                      interaction_matrix::DefaultDict{Tuple{String, String}, Bool, Bool}, box::Vector{Float64}, multithreaded::Bool)
    ϵ = Dict{Tuple{String, String}, Float64}()
    σ = Dict{Tuple{String, String}, Float64}()
    r_cut = Dict{Tuple{String, String}, Float64}()
    for (id1, id2, vals) in params
        ϵ[id1, id2] = ϵ[id2, id1] = vals["ϵ"]
        σ[id1, id2] = σ[id2, id1] = vals["σ"]
        r_cut[id1, id2] = r_cut[id2, id1] = vals["r_cut"]
    end

    return LennardJones(ϵ, σ, r_cut, particles, neighbor_list, interaction_matrix, box, multithreaded)
end

function compute_forces!(lj::LennardJones)
    @use_threads lj.multithreaded for particle in lj.particles
        x, y, z = particle.position
        i = floor(Int64, mod(x, lj.box[1]) / lj.neighbor_list.cell_spacings[1])
        j = floor(Int64, mod(y, lj.box[2]) / lj.neighbor_list.cell_spacings[2])
        k = floor(Int64, mod(z, lj.box[3]) / lj.neighbor_list.cell_spacings[3])
        for Δi = -1 : 1, Δj = -1 : 1
            iΔi = mod(i + Δi, lj.neighbor_list.cell_counts[1]) + 1
            jΔj = mod(j + Δj, lj.neighbor_list.cell_counts[2]) + 1
            kΔk = mod(k + Δk, lj.neighbor_list.cell_counts[3]) + 1
            index = lj.neighbor_list.start_index[iΔi, jΔj, kΔk]
            while id > 0
                neighbor = lj.neighbor_list.particles[index]
                if isnothing(particle.body_id) || isnothing(neighbor.body_id) || particle.body_id != neighbor.body_id
                    Δx, Δy, Δz = wrap_displacement(x - neighbor.position[1], y - neighbor.position[2], z - neighbor.position[3], lj.box)
                    
                    Δr² = Δx^2 + Δy^2 + Δz^2
                    if lj.interaction_matrix[particle.id, neighbor.id] && 0.0 < Δr² < lj.r_cut[particle.id, neighbor.id]^2
                        val = (lj.σ[particle.id, neighbor.id]^2 / Δr²)^3
                        coef = lj.ϵ[particle.id, neighbor.id] * (48.0 * val - 24.0) * val / Δr²

                        particle.force[1] += coef * Δx
                        particle.force[2] += coef * Δy
                        particle.force[3] += coef * Δz
                    end
                    index = lj.neighbor_list.next_index[index]
                end
            end
        end
    end
end

###########################################################################################################################################
# 
###########################################################################################################################################

mutable struct Morse <: Interaction
    D0::Dict{Tuple{String, String}, Float64}
    α::Dict{Tuple{String, String}, Float64}
    r0::Dict{Tuple{String, String}, Float64}
    r_cut::Dict{Tuple{String, String}, Float64}

    particles::Vector{Particle}
    neighbor_list::LinkedCellList
    interaction_matrix::DefaultDict{Tuple{String, String}, Bool, Bool}
    box::SVector{3, Float64}

    multithreaded::Bool
end

function Morse(params::Vector{Tuple{String, String, Dict{String, Float64}}}, particles::Vector{Particle}, neighbor_list::LinkedCellList, 
               interaction_matrix::DefaultDict{Tuple{String, String}, Bool, Bool}, box::Vector{Float64}, multithreaded::Bool)
    D0 = Dict{Tuple{String, String}, Float64}()
    α = Dict{Tuple{String, String}, Float64}()
    r0 = Dict{Tuple{String, String}, Float64}()
    r_cut = Dict{Tuple{String, String}, Float64}()
    for (id1, id2, vals) in params
        D0[id1, id2] = D0[id2, id1] = vals["D0"]
        α[id1, id2] = α[id2, id1] = vals["α"]
        r0[id1, id2] = r0[id2, id1] = vals["r0"]
        r_cut[id1, id2] = r_cut[id2, id1] = vals["r_cut"]
    end

    return Morse(D0, α, r0, r_cut, particles, neighbor_list, interaction_matrix, box, multithreaded)
end

function compute_forces!(m::Morse)
    @use_threads m.multithreaded for particle in m.particles
        x, y, z = particle.position
        i = floor(Int64, mod(x, m.box[1]) / m.neighbor_list.cell_spacings[1])
        j = floor(Int64, mod(y, m.box[2]) / m.neighbor_list.cell_spacings[2])
        k = floor(Int64, mod(z, m.box[3]) / m.neighbor_list.cell_spacings[3])
        for Δi = -1 : 1, Δj = -1 : 1
            iΔi = mod(i + Δi, m.neighbor_list.cell_counts[1]) + 1
            jΔj = mod(j + Δj, m.neighbor_list.cell_counts[2]) + 1
            kΔk = mod(k + Δk, m.neighbor_list.cell_counts[3]) + 1
            index = m.neighbor_list.start_index[iΔi, jΔj, kΔk]
            while id > 0
                neighbor = m.neighbor_list.particles[index]
                if isnothing(particle.body_id) || isnothing(neighbor.body_id) || particle.body_id != neighbor.body_id
                    Δx, Δy, Δz = wrap_displacement(x - neighbor.position[1], y - neighbor.position[2], z - neighbor.position[3], m.box)
                    
                    Δr² = Δx^2 + Δy^2 + Δz^2
                    if m.interaction_matrix[particle.id, neighbor.id] && 0.0 < Δr² < m.r_cut[particle.id, neighbor.id]^2
                        Δr = sqrt(Δr²)
                        val = exp(-m.α[particle.id, neighbor.id] * (Δr - m.r0[particle.id, neighbor.id]))
                        coef = 2 * m.α[particle.id, neighbor.id] * m.D0[particle.id, neighbor.id] * val * (val - 1.0) / Δr

                        particle.force[1] += coef * Δx
                        particle.force[2] += coef * Δy
                        particle.force[3] += coef * Δz
                    end
                    index = m.neighbor_list.next_index[index]
                end
            end
        end
    end
end

###########################################################################################################################################
# 
###########################################################################################################################################

mutable struct HarmonicBond <: Interaction
    k::Float64
    r0::SVector{3, Float64}

    bond_list::BondList
    box::SVector{3, Float64}

    multithreaded::Bool
end

function compute_forces!(hb::HarmonicBond)
    @use_threads hb.multithreaded for (particle, neighbors) in hb.bond_list
        for neighbor in neighbors
            Δx, Δy, Δz = wrap_displacement(particle.position[1] - neighbor.position[1], particle.position[2] - neighbor.position[2], particle.position[3] - neighbor.position[3], hb.box)
            
            coef = -hb.k * (hb.r0 == 0.0 ? 1.0 : 1.0 - hb.r0 / sqrt(Δx^2 + Δy^2 + Δz^2))
            particle.force[1] += coef * Δx
            particle.force[2] += coef * Δy
            particle.force[3] += coef * Δz
        end
    end
end

###########################################################################################################################################
# Additional functions
###########################################################################################################################################

function wrap_displacement(Δx::Float64, Δy::Float64, Δz::Float64, box::SVector{3, Float64})
    if abs(Δx) > 0.5 * box[1]
        Δx -= sign(Δx) * box[1]
    end
    if abs(Δy) > 0.5 * box[2]
        Δy -= sign(Δy) * box[2]
    end
    if abs(Δz) > 0.5 * box[3]
        Δz -= sign(Δz) * box[3]
    end

    return Δx, Δy, Δz
end

function create_interaction_matrix(pairs::Vector{Tuple{String, String}})
    interaction_matrix = DefaultDict{Tuple{String, String}, Bool}(false)
    for (id1, id2) in pairs
        interaction_matrix[(id1, id2)] = interaction_matrix[(id2, id1)] = true
    end

    return interaction_matrix
end