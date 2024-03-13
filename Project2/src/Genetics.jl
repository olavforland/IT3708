module Genetics

export Chromosome, compute_fitness!, compute_unfitness!, simple_routing_sort_start_times

using ..DataParser: Patient, ProblemInstance

mutable struct Chromosome

    genotype::Vector{Int}
    phenotype::Union{Vector{Vector{Int}},Nothing}

    fitness::Union{Float64,Nothing}
    time_unfitness::Union{Float64,Nothing}
    strain_unfitness::Union{Float64,Nothing}

    route_fitness::Union{Vector{Float64},Nothing}
    route_time_unfitness::Union{Vector{Float64},Nothing}
    route_strain_unfitness::Union{Vector{Float64},Nothing}

    # Constructor that sets phenotype, fitness and unfitness to nothing by default
    Chromosome(genotype::Vector{Int}) = new(genotype, nothing, nothing, nothing, nothing, nothing, nothing)
    Chromosome(genotype::Vector{Int}, n_nurses::Int) = new(
        genotype,
        [Vector{Int}() for _ in 1:n_nurses],
        nothing, nothing, nothing,
        [0.0 for _ in 1:n_nurses],
        [0.0 for _ in 1:n_nurses],
        [0.0 for _ in 1:n_nurses]
    )

    id::Int
    patient_sets::Vector{Set{Int}}

end




function compute_fitness!(chromosome::Chromosome, problem_instance::ProblemInstance)
    travel_times = problem_instance.travel_times
    travel_time = 0.0
    for (n, route) in enumerate(chromosome.phenotype)
        route_travel_time = 0.0
        # Start at depot
        prev_patient = 1
        for patient in route
            # Increment by one due to 1-indexing
            route_travel_time += travel_times[prev_patient][patient+1]
            # Increment by one due to 1-indexing
            prev_patient = patient + 1
        end
        # End at depot
        route_travel_time += travel_times[prev_patient][1]
        chromosome.route_fitness[n] = route_travel_time
        travel_time += route_travel_time

    end

    chromosome.fitness = travel_time
end

function compute_unfitness!(chromosome::Chromosome, problem_instance::ProblemInstance)
    time_unfitness = 0.0
    strain_unfitness = 0.0

    for (n, route) in enumerate(chromosome.phenotype)
        elapsed_time = 0.0
        time_violation = 0.0
        nurse_strain = 0
        prev_patient = 1
        for patient in route
            # Add travel time
            elapsed_time += problem_instance.travel_times[prev_patient][patient+1]
            # If arrive early, wait
            elapsed_time += max(problem_instance.patients[patient].start_time - elapsed_time, 0)
            # If arrive late, add to time violation
            time_violation += max(elapsed_time - (problem_instance.patients[patient].end_time - problem_instance.patients[patient].care_time), 0)
            # Add care time
            elapsed_time += problem_instance.patients[patient].care_time
            # Add strain on nurse
            nurse_strain += problem_instance.patients[patient].demand

            prev_patient = patient + 1
        end
        # Add travel time back to depot
        elapsed_time += problem_instance.travel_times[prev_patient][1]
        time_violation += max(elapsed_time - problem_instance.depot_return_time, 0)

        # Express time unfitness as a percentage of the total time spent
        route_time_unfitness = time_violation > 0 ? time_violation / elapsed_time : 0.0
        chromosome.route_time_unfitness[n] = route_time_unfitness
        time_unfitness += route_time_unfitness

        # Express strain unfitness as a percentage of the total capacity
        route_strain_unfitness = max(nurse_strain - problem_instance.nurse_capacity, 0) / problem_instance.nurse_capacity
        chromosome.route_strain_unfitness[n] = route_strain_unfitness
        strain_unfitness += route_strain_unfitness
    end

    chromosome.time_unfitness = time_unfitness
    chromosome.strain_unfitness = strain_unfitness
end


function simple_routing_sort_start_times(chromosome::Chromosome, patients::Vector{Patient}, n_nurses::Int)
    routes = [Vector{Int}() for _ in 1:n_nurses]  # Initialize empty routes for each nurse

    # Populate routes with patient indices
    for (patient_index, nurse_id) in enumerate(chromosome.genotype)
        push!(routes[nurse_id], patient_index)
    end

    # Sort routes based on patients' start times
    for i in 1:n_nurses
        routes[i] = sort(routes[i], by=p -> patients[p].start_time)
    end

    return routes
end

end # module