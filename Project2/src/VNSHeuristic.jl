module VNSHeuristic

export construct_solution!, improve_solution!

using Random
using ..DataParser: ProblemInstance, Patient
using ..Genetics: Chromosome

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
    solution, obj = local_1_shift(solution, instance, construction_objective)
    # Iterate until feasibility or max level is reached
    while obj > 0 && level < 5
        # Perform level random 1-shifts
        new_solution = random_1_shift(solution, level)
        # Perform local search
        new_solution, new_obj = local_1_shift(new_solution, instance, construction_objective)
        # Save new solution if better
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
    route = variable_neighborhood_decent(route, instance, improvement_objective)
    while (level <= 5) && (iter < max_iter)

        new_route = random_1_shift(route, level)    
        new_route = variable_neighborhood_decent(new_route, instance, improvement_objective)

        if improvement_objective(new_route, instance) < improvement_objective(route, instance)
            level = 1
            route = new_route
        else
            level += 1
        end
        iter += 1
    end


    return route

end

function variable_neighborhood_decent(route::Vector{Patient}, instance::ProblemInstance, objective::Function)
    prev_route = nothing

    max_iter = 10
    iter = 0
    
    while (prev_route != map(p -> p.id, route)) && (iter < max_iter)
        route, _ = local_1_shift(route, instance, objective)
        prev_route = map(p -> p.id, route)
        
        route = local_2_opt(route, instance, objective)

        iter += 1

    end
    return route

end

function local_2_opt(route::Vector{Patient}, instance::ProblemInstance, objective::Function)
    best_route = route
    best_obj = objective(route, instance)

    for i in 1:length(route)
        for j in i+1:length(route)
            # If j cannot precede i, skip
            if (route[j].id, route[i].id) ∈ instance.inadmissable_presedence
                break
            end
            new_route = deepcopy(route)

            new_route[i:j] = reverse(new_route[i:j])
            
            obj = objective(new_route, instance)
            if obj < best_obj
                best_obj = obj
                best_route = new_route
            end
        end
    end

    return best_route

end

# ---------------- Helpers ----------------

function local_1_shift(route::Vector{Patient}, instance::ProblemInstance, objective::Function)
    
    best_route = route
    best_obj = objective(route, instance)
    
    # Iterate through all possible 1-shifts
    for i in 1:length(route)
        # Consider the cases when patient i precedes patient j
        for j in i+1:length(route)
            # If j cannot precede i, i cannot be inserted after j
            if (route[j].id, route[i].id) ∈ instance.inadmissable_presedence
                break
            end

            shifted_patient = splice!(route, i) # Remove the element at position j
            insert!(route, j - 1, shifted_patient) # Insert that element at position i
            # Recalculate objective
            obj = objective(route, instance)
            if obj < best_obj
                best_obj = obj
                best_route = route
            end
        end
        # Consider the cases when patient i succeeds patient j
        for j in i-1:-1:1
            # If j cannot precede i, i cannot be inserted after j
            if (route[i].id, route[j].id) ∈ instance.inadmissable_presedence
                break
            end

            shifted_patient = splice!(route, i) # Remove the element at position j
            insert!(route, j, shifted_patient) # Insert that element at position i
            # Recalculate objective
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
    # Perform level random 1-shifts 
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

function improvement_objective(route::Vector{Patient}, instance::ProblemInstance)

    travel_times = instance.travel_times
    travel_time = 0.0
    # Start at depot
    prev_patient = 1
    for patient in route
        # Increment by one due to 1-indexing
        travel_time += travel_times[prev_patient][patient.id + 1]
        # Increment by one due to 1-indexing
        prev_patient = patient.id + 1
    end
    # End at depot
    travel_time += travel_times[prev_patient][1]

    return travel_time

end

end # module