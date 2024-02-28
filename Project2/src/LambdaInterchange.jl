module LambdaInterchange

export shift_operation

using ..Genetics: Chromosome
using ..Objective: time_unfitness_objective, time_fitness_objective, strain_unfitness_objective
using ..Genetics: compute_fitness!, compute_unfitness!
using ..DataParser: ProblemInstance, Patient

function shift_operation(chromosome::Chromosome, problem_instance::ProblemInstance)


    improved = true

    best_chromosome = deepcopy(chromosome)
    iter = 0
    while improved && iter < 100
        improved = false
        iter += 1
        for i in 1:length(chromosome.phenotype)
            for j in 1:length(chromosome.phenotype)
                if i == j # Ensure we're not trying to shift within the same route
                    continue
                end
                for patient in best_chromosome.phenotype[i]
                    new_chromosome = deepcopy(best_chromosome)

                    deleteat!(new_chromosome.phenotype[i], findfirst(==(patient), new_chromosome.phenotype[i]))
                    push!(new_chromosome.phenotype[j], patient)

                    compute_fitness!(new_chromosome, problem_instance)
                    compute_unfitness!(new_chromosome, problem_instance)

                    if (new_chromosome.fitness <= best_chromosome.fitness) && 
                       (new_chromosome.time_unfitness <= best_chromosome.time_unfitness) && 
                       (new_chromosome.strain_unfitness <= best_chromosome.strain_unfitness)

                        best_chromosome = new_chromosome
                        improved = true
                        break  # Exit the innermost loop after improvement
                    end
                end
                if improved
                    break  # Exit the middle loop after improvement
                end
            end
            if improved
                break  # Exit the outermost loop after improvement
            end
        end
    end
    return best_chromosome
end


end # module