module Mutation

export swap_mutation!, insert_mutation!

using ..Genetics: Chromosome



function swap_mutation!(chromosome::Chromosome, mutation_rate::Float64=0.1)

    if rand() > mutation_rate
        return
    end

    n_patients = length(chromosome.genotype)

    # Randomly select two patients
    patient1 = rand(1:n_patients)
    patient2 = rand(1:n_patients)

    # Swap the patients
    chromosome.genotype[patient1], chromosome.genotype[patient2] = chromosome.genotype[patient2], chromosome.genotype[patient1]

    # Phenotype, fitness and unfitness must be recomputed
    chromosome.phenotype = [Vector{Int}() for _ in 1:length(chromosome.phenotype)]
    chromosome.fitness = nothing
    chromosome.time_unfitness = nothing
    chromosome.strain_unfitness = nothing
end

function insert_mutation!(chromosome::Chromosome, mutation_rate::Float64=0.1)

    if rand() > mutation_rate
        return
    end

    n_patients = length(chromosome.genotype)

    # Randomly select two patients
    patient1 = rand(1:n_patients)
    patient2 = rand(1:n_patients)

    n_used_patients = length(unique(chromosome.genotype))

    # Find a new nurse and assign to it patient i with probability 1/n_used_patients
    new_nurse = false
    for i in minimum(chromosome.genotype):maximum(chromosome.genotype)
        if i ∉ chromosome.genotype
            if rand() > 1 / n_used_patients
                break
            end

            chromosome.genotype[patient1] = i
            new_nurse = true
            break
        end
    end

    # Assign patient 1 to nurse handling patient 2 if not new nurse is used
    if !new_nurse
        chromosome.genotype[patient1] = chromosome.genotype[patient2]
    end

    # Phenotype, fitness and unfitness must be recomputed
    chromosome.phenotype = [Vector{Int}() for _ in 1:length(chromosome.phenotype)]
    chromosome.fitness = nothing
    chromosome.time_unfitness = nothing
    chromosome.strain_unfitness = nothing
end

end # module