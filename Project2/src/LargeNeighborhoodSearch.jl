module LargeNeighborhoodSearch

using ..DataParser: ProblemInstance, Patient
using ..Genetics: Chromosome, compute_fitness!, compute_unfitness!
using ..LocalSearch
using ..Objective


# Solves a TSP for all routes in the chromosome, with the depot replicated between routes
function tsp_all_routes!(chromosome::Chromosome, problem_instance::ProblemInstance)

    # println("Starting tsp_all_routes")

    complete_route = Vector{Int}()


    for (i, route) in enumerate(chromosome.phenotype)
        # Append route to result
        append!(complete_route, route)

        # Add depot between each vector, but not after last one
        if i < length(chromosome.phenotype)
            push!(complete_route, 0)
        end
    end
        
    # Convert the route to patients, with special care for depot
    patients = Vector{Patient}()
    nurse_nr = 1
    return_time = problem_instance.depot_return_time
    for patient_id in complete_route
        if patient_id == 0
            x, y = problem_instance.depot_coords
            push!(patients, Patient(0, x, y, 0, 0, problem_instance.depot_return_time * problem_instance.n_nurses, 0))
            nurse_nr += 1
        else
            original_patient = problem_instance.patients[patient_id]
            # Push each time window depot_return_time * nurse_nr to ensure feasible solutions
            augmented_patient = Patient(
                patient_id, 
                original_patient.x_coord, 
                original_patient.y_coord, 
                original_patient.demand, 
                original_patient.start_time + return_time * (nurse_nr - 1), 
                original_patient.end_time + return_time * (nurse_nr - 1), 
                original_patient.care_time, 
            )
            push!(patients, augmented_patient)
        end
    end
        
    local_2_opt!(patients, problem_instance, total_objective)

    routes = Vector{Vector{Int}}()
    prev_zero_ind = 1
    for (i, el) in enumerate(patients)
        if el.id == 0
            push!(routes, map(p -> p.id, patients[prev_zero_ind:i-1]))
            prev_zero_ind = i + 1
        end
    end
    # Push last route
    push!(routes, map(p -> p.id, patients[prev_zero_ind:end]))
    
    while length(routes) < problem_instance.n_nurses
        push!(routes, Vector{Int}())  # Add an empty Vector{Int}
    end


    chromosome.phenotype = routes
    for (nurse, route) in enumerate(routes)
        for patient in route
            chromosome.genotype[patient] = nurse
        end
    end

    compute_fitness!(chromosome, problem_instance)
    compute_unfitness!(chromosome, problem_instance)

    # println("Finished tsp_all_routes")

end

end # module