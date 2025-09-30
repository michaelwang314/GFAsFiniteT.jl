struct System
    box::SVector{3, Float64}

    bodies::Vector{<:Body}
    interactions::Vector{<:Interaction}
    integrator::Integrator
end