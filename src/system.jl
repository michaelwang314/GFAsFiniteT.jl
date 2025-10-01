struct System
    box::SVector{3, Float64}

    bodies::Vector{<:Body}
    interactions::Vector{<:Interaction}
    external_forces::Vector{<:ExternalForce}
    integrator::Integrator
end