module LambdaInterchange

export lambda_shift_operation, lambda_interchange_operation

using ..Genetics: Chromosome
using ..Objective: time_unfitness_objective, time_fitness_objective, strain_unfitness_objective, total_objective
using ..Genetics: compute_fitness!, compute_unfitness!
using ..DataParser: ProblemInstance, Patient
using ..VNSHeuristic: construct_single_route, local_2_opt!

function lambda_shift_operation(chromosome::Chromosome, problem_instance::ProblemInstance)

    improved = true

    best_chromosome = deepcopy(chromosome)

    while improved
        improved = false

        # None empty indices of chromosome.phenotype
        non_empty_indices = [index for (index, vector) in enumerate(chromosome.phenotype) if !isempty(vector)]

        # Find the index of the first empty vector
        first_empty_index = findfirst(isempty, chromosome.phenotype)
        # Check if a first empty vector was found and append its index if so
        nurse_indices = first_empty_index !== nothing ? [non_empty_indices; first_empty_index] : non_empty_indices

        for i in nurse_indices
            for j in nurse_indices
                if i == j # Ensure we're not trying to shift within the same route
                    continue
                end
                for patient in best_chromosome.phenotype[i]

                    # Save previous fitness and unfitness values
                    best_fitness = best_chromosome.route_fitness[i] + best_chromosome.route_fitness[j]
                    best_time_unfitness = best_chromosome.route_time_unfitness[i] + best_chromosome.route_time_unfitness[j]
                    best_strain_unfitness = best_chromosome.route_strain_unfitness[i] + best_chromosome.route_strain_unfitness[j]
                    
                    # Assign patient to new nurse
                    best_chromosome.genotype[patient] = j

                    improved = evaluate_lambda_improvement!(best_chromosome, i, j, problem_instance, best_fitness, best_time_unfitness, best_strain_unfitness)

                    if improved
                        break  # Exit the innermost loop after improvement
                    else
                        best_chromosome.genotype[patient] = i
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

function lambda_interchange_operation(chromosome::Chromosome, problem_instance::ProblemInstance)
    improved = true

    best_chromosome = deepcopy(chromosome)

    while improved
        improved = false
        # None empty indices of chromosome.phenotype
        non_empty_indices = [index for (index, vector) in enumerate(best_chromosome.phenotype) if !isempty(vector)]

        # Find the index of the first empty vector
        first_empty_index = findfirst(isempty, best_chromosome.phenotype)
        # Check if a first empty vector was found and append its index if so
        nurse_indices = first_empty_index !== nothing ? [non_empty_indices; first_empty_index] : non_empty_indices

        for i in nurse_indices
            for j in nurse_indices
                if i == j # Skip if it's the same route
                    continue
                end

                for patient_i in best_chromosome.phenotype[i]
                    for patient_j in best_chromosome.phenotype[j]

                        # Save previous fitness and unfitness values
                        best_fitness = best_chromosome.route_fitness[i] + best_chromosome.route_fitness[j]
                        best_time_unfitness = best_chromosome.route_time_unfitness[i] + best_chromosome.route_time_unfitness[j]
                        best_strain_unfitness = best_chromosome.route_strain_unfitness[i] + best_chromosome.route_strain_unfitness[j]

                        # Swap the patients between the two nurses/routes
                        best_chromosome.genotype[patient_i], best_chromosome.genotype[patient_j] = best_chromosome.genotype[patient_j], best_chromosome.genotype[patient_i]

                        improved = evaluate_lambda_improvement!(best_chromosome, i, j, problem_instance, best_fitness, best_time_unfitness, best_strain_unfitness)

                        if improved
                            break  # Exit the innermost loop after improvement
                        else
                            # Swap the patients back
                            best_chromosome.genotype[patient_i], best_chromosome.genotype[patient_j] = best_chromosome.genotype[patient_j], best_chromosome.genotype[patient_i]
                        end
                    end
                    if improved
                        break  # Break from the second loop after improvement
                    end
                end
                if improved
                    break  # Break from the third loop after improvement
                end
            end
        end
    end

    return best_chromosome
end

# ---------------- Helpers ----------------

function evaluate_lambda_improvement!(chromosome::Chromosome, i::Int, j::Int, problem_instance::ProblemInstance, prev_fitness::Float64, prev_time_unfitness::Float64, prev_strain_unfitness::Float64)

    nurse_patients_i = findall(x -> x == i, chromosome.genotype)
    nurse_patients_j = findall(x -> x == j, chromosome.genotype)

    route_i = construct_single_route(problem_instance, problem_instance.patients[nurse_patients_i])
    route_j = construct_single_route(problem_instance, problem_instance.patients[nurse_patients_j])

    time_unfitness_i, strain_unfitness_i, fitness_i = total_objective(route_i, problem_instance)
    time_unfitness_j, strain_unfitness_j, fitness_j = total_objective(route_j, problem_instance)

    if  (fitness_i + fitness_j <= prev_fitness) && 
        (time_unfitness_i + time_unfitness_j <= prev_time_unfitness) &&
        (strain_unfitness_i + strain_unfitness_j <= prev_strain_unfitness) &&
        (fitness_i + fitness_j < prev_fitness || time_unfitness_i + time_unfitness_j < prev_time_unfitness || strain_unfitness_i + strain_unfitness_j < prev_strain_unfitness)

        # Improve promising routes with 2-opt
        local_2_opt!(route_i, problem_instance, total_objective)
        local_2_opt!(route_j, problem_instance, total_objective)

        # Add routes to best chromosome
        chromosome.phenotype[i] = map(p -> p.id, route_i)
        chromosome.phenotype[j] = map(p -> p.id, route_j)

        # Update fitness and unfitness values
        chromosome.route_fitness[i] = fitness_i
        chromosome.route_fitness[j] = fitness_j

        chromosome.route_time_unfitness[i] = time_unfitness_i
        chromosome.route_time_unfitness[j] = time_unfitness_j

        chromosome.route_strain_unfitness[i] = strain_unfitness_i
        chromosome.route_strain_unfitness[j] = strain_unfitness_j

        chromosome.fitness += - prev_fitness + fitness_i + fitness_j
        chromosome.time_unfitness += - prev_time_unfitness + time_unfitness_i + time_unfitness_j
        chromosome.strain_unfitness += - prev_strain_unfitness + strain_unfitness_i + strain_unfitness_j
        
        # Sometimes numerical instable
        chromosome.time_unfitness = max(chromosome.time_unfitness, 0)
        chromosome.strain_unfitness = max(chromosome.strain_unfitness, 0)

        return true
    end
    return false

end

    

end # module