module Genetics

export Chromosome, compute_fitness!

using ..DataParser: Patient, ProblemInstance

mutable struct Chromosome
    
    genotype::Vector{Int}
    phenotype::Union{Vector{Vector{Int}}, Nothing}
    
    fitness::Union{Float64, Nothing}
    unfitness::Union{Float64, Nothing}

    # Constructor that sets phenotype, fitness and unfitness to nothing by default
    Chromosome(genotype::Vector{Int}) = new(genotype, nothing, nothing, nothing)
end


function compute_fitness!(chromosome::Chromosome, problem_instance::ProblemInstance)
    travel_times = problem_instance.travel_times
    travel_time = 0.0
    for route in chromosome.phenotype

        # Start at depot
        prev_patient = 1
        for patient in route
            travel_time += travel_times[prev_patient][patient]
            prev_patient = patient
        end
        # End at depot
        travel_time += travel_times[prev_patient][1]
    end

    chromosome.fitness = travel_time
end

end # module