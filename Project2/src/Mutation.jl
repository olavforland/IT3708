module Mutation

export swap_mutation!

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

end # module