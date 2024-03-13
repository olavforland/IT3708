module VNSHeuristic

export construct_solution!, improve_solution!, construct_single_route, improve_single_route

using Random
using ..DataParser: ProblemInstance, Patient
using ..Genetics: Chromosome
using ..Objective: time_unfitness_objective, time_fitness_objective, total_objective
using ..LocalSearch

# ---------------- Construction Heuristic ----------------

function construct_solution!(instance::ProblemInstance, chromosome::Chromosome, n_nurses::Int)
    for nurse in 1:n_nurses
        nurse_patients = findall(x -> x == nurse, chromosome.genotype)
        route = construct_single_route(instance, instance.patients[nurse_patients])
        chromosome.phenotype[nurse] = map(p -> p.id, route)
    end
end


function construct_single_route(instance::ProblemInstance, patients::Vector{Patient})

    level = 1
    
    # Start from a random initial solution
    solution = shuffle(patients)
    # Perform local search
    obj = local_1_shift!(solution, instance, time_unfitness_objective)

    # Iterate until feasibility or max level is reached
    while obj > 0 && level < 5
        # Perform level random 1-shifts
        swaps = random_1_shift!(solution, level)
        # Perform local search
        new_obj = local_1_shift!(solution, instance, time_unfitness_objective)
        # Save new solution if better
        if new_obj < obj
            level = 1
            obj = new_obj
        else
            level += 1
            # Swap back to original positions
            for swap in reverse(swaps)
                i, j = swap
                solution[i], solution[j] = solution[j], solution[i] 
            end
        end
    end

    return solution
end

# ---------------- Improvement Heuristic ----------------

function improve_solution!(instance::ProblemInstance, chromosome::Chromosome)
    for nurse in 1:instance.n_nurses
        nurse_patients = findall(x -> x == nurse, chromosome.genotype)
        # Only improve routes with more than 3 patients
        if length(nurse_patients) <= 4
            continue
        end

        route = improve_single_route(instance.patients[nurse_patients], instance)
        chromosome.phenotype[nurse] = map(p -> p.id, route)
        
    end
end


function improve_single_route(route::Vector{Patient}, instance::ProblemInstance)
    level = 1
    max_iter = 10
    iter = 0
    best_route = deepcopy(route)
    best_obj = variable_neighborhood_decent!(best_route, instance, total_objective)
    while (level <= 5) && (iter < max_iter)

        new_route = deepcopy(best_route)

        random_1_shift!(new_route, level)
        
        obj = variable_neighborhood_decent!(new_route, instance, total_objective)

        if obj < best_obj
            level = 1
            best_route = deepcopy(new_route)
        else
            level += 1
        end
        iter += 1
    end


    return best_route

end


end # module