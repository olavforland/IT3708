module VNSHeuristic

export construct_solution!

using Random
using ..DataParser: ProblemInstance, Patient
using ..Genetics: Chromosome


function construct_solution!(instance::ProblemInstance, chromosome::Chromosome)
    for nurse in 1:instance.n_nurses
        nurse_patients = findall(x -> x == nurse, chromosome.genotype)
        route = construct_single_route(instance, instance.patients[nurse_patients])
        chromosome.phenotype[nurse] = map(p -> p.id, route)
    end
end


function construct_single_route(instance::ProblemInstance, patients::Vector{Patient})

    level = 1
    
    solution = shuffle(patients)#sort(patients, by = p -> p.start_time)
    solution, obj = local_1_shift(solution, instance, construction_objective)

    while obj > 0 && level < 10
        new_solution = random_1_shift(solution, level)
        new_solution, new_obj = local_1_shift(new_solution, instance, construction_objective)
        
        if new_obj < obj
            level = 1
            solution = new_solution
            obj = new_obj
        else
            level += 1
        end
    end

    return solution
end


# ---------------- Helpers ----------------

function local_1_shift(route::Vector{Patient}, instance::ProblemInstance, objective::Function)
    
    best_route = route
    best_obj = objective(route, instance)
    
    
    for i in 1:length(route)
        for j in i+1:length(route)
            if (route[j].id, route[i].id) ∈ instance.inadmissable_presedence
                break
            end
            route[i], route[j] = route[j], route[i]
            obj = objective(route, instance)
            if obj < best_obj
                best_obj = obj
                best_route = route
            end
        end
    end
    
    return best_route, best_obj
end

function random_1_shift(route::Vector{Patient}, level::Int)
    for _ in 1:level
        i = rand(1:length(route))
        j = rand(1:length(route))
        if i != j
            route[i], route[j] = route[j], route[i]
        end
    end
    
    return route
    
end

function construction_objective(route::Vector{Patient}, instance::ProblemInstance)

    elapsed_time = 0.0
    time_violation = 0.0

    prev_patient = 1
    for patient in route
        # Add travel time
        elapsed_time += instance.travel_times[prev_patient][patient.id + 1]
        # If arrive early, wait
        elapsed_time += max(patient.start_time - elapsed_time, 0)
        # If arrive late, add to time violation
        time_violation += max(elapsed_time - (patient.end_time - patient.care_time), 0)
        # Add care time
        elapsed_time += patient.care_time
        # Increment by one due to 1-indexing
        prev_patient = patient.id + 1
    end

    return time_violation
end

end # module