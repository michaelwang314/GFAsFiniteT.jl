module GFAsFiniteT
    using StaticArrays
    using Serialization
    using LinearAlgebra
    using DataStructures

    export Body, Particle, RigidBody, rotate!, translate!, set_body_ids!, rigid_bodies_to_particle_list, get_particles_with_ids
    export NeighborList, BondList, LinkedCellList, bind_closest, update_cell_list!
    export Interaction, LennardJones, Morse, HarmonicBond, compute_forces!, wrap_displacement, create_interaction_matrix
    export ExternalForce, ConstantForce
    export Integrator, Brownian, update_bodies!
    export System, Trajectories, run_simulation!, export_for_mathematica!
    export use_threads

    """
    Flush output so that jobs can be monitored on a cluster.
    """
    @inline println(args...) = println(stdout, args...)
    @inline function println(io::IO, args...)
        Base.println(io, args...)
        flush(io)
    end

    macro use_threads(multithreaded::Union{Expr, Symbol}, expr::Expr)
        esc(quote
            if $multithreaded
                Threads.@threads $expr
            else
                $expr
            end
        end)
    end

    include("bodies.jl")
    include("neighborlists.jl")
    include("interactions.jl")
    include("externalforces.jl")
    include("integrators.jl")
    include("system.jl")
end
