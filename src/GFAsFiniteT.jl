module GFAsFiniteT
    using StaticArrays
    using Serialization
    using LinearAlgebra
    using DataStructures

    export Body, Particle, RigidBody, rotate!, translate!, set_body_ids!, get_particle_list, get_particles_with_ids
    export NeighborList, BondList, LinkedCellList, bind_closest, update_neighbor_list!
    export Interaction, LennardJones, Morse, HarmonicBond, compute_forces!, correct_for_periodicity, create_interaction_matrix
    export ExternalForce, ConstantForce, HarmonicTrap
    export Integrator, Brownian, update_bodies!
    export System, Trajectories, run_simulation!, initialize_random_positions, save!, load
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

    #include("visualizer/Visualizer.jl")
end
