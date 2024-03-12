module VismaCrossover


using ..DataParser: ProblemInstance
using ..VNSHeuristic: construct_solution!, improve_solution!
using ..Genetics: Chromosome
using ..Objective: total_objective


export visma_crossover



function visma_crossover(p1::Chromosome, p2::Chromosome, n_nurses::Int, instance::ProblemInstance)
    """
    1: Select a nurse from both parents.
    2: Remove the nurses patients from the other parent.
    3: For each patient without a nurse, try to assign it to all other nurses and reoptimize the routes.
    4: Assign the patient to the nurse with the best route, after (time_unfitness, strain_unfitness, fitness)
    """
    n_patients = length(p1.genotype)

    # Create children
    c1 = Chromosome([0 for _ in 1:n_patients], n_nurses)
    c2 = Chromosome([0 for _ in 1:n_patients], n_nurses)

    # Randomly select two nurses
    nurses_p1 = unique(p1.genotype)
    nurses_p2 = unique(p2.genotype)

    nurse_p1 = nurses_p1[rand(1:length(nurses_p1))]
    nurse_p2 = nurses_p2[rand(1:length(nurses_p2))]

    #find the patients for the selected nurses
    nurse_patients_p1 = findall(x -> x == nurse_p1, p1.genotype)

    nurse_patients_p2 = findall(x -> x == nurse_p2, p2.genotype)


    #Find one unused nurse in either parent if there are any
    all_nurses = Set(1:n_nurses)
    unused_nurses_p1 = setdiff(all_nurses, unique(p1.genotype))
    unused_nurses_p2 = setdiff(all_nurses, unique(p2.genotype))


    for i in 1:n_patients
        c1.genotype[i] = p1.genotype[i]
        c2.genotype[i] = p2.genotype[i]
    end

    #Add one unused nurse to the nurses in parent
    if !isempty(unused_nurses_p1)
        push!(nurses_p1, first(unused_nurses_p1))
    end
    if !isempty(unused_nurses_p2)
        push!(nurses_p2, first(unused_nurses_p2))
    end

    # Try for all patients in nurse_patients_p1, try to assign them to a different nurse in c1
    reallocations_c1 = Vector{Int}()
    for patient in nurse_patients_p2
        best_nurse = 0
        best_obj = (Inf, Inf, Inf)
        curr_nurse = c1.genotype[patient]
        for nurse in nurses_p1
            if nurse != curr_nurse
                c1.genotype[patient] = nurse
                construct_solution!(instance, c1, n_nurses)
                route = [instance.patients[p] for p in c1.phenotype[nurse]]
                obj = total_objective(route, instance)
                if obj[1] < best_obj[1] ||
                   (obj[1] == best_obj[1] && obj[2] < best_obj[2]) ||
                   (obj[1] == best_obj[1] && obj[2] == best_obj[2] && obj[3] < best_obj[3])
                    best_nurse = nurse
                    best_obj = obj
                end
            end
        end
        c1.genotype[patient] = best_nurse
        reallocations_c1 = push!(reallocations_c1, best_nurse)
    end

    # Try for all patients in nurse_patients_p2, try to assign them to a different nurse in c2
    reallocations_c2 = Vector{Int}()
    for patient in nurse_patients_p1
        best_nurse = 0
        best_obj = (Inf, Inf, Inf)
        curr_nurse = c2.genotype[patient]
        for nurse in nurses_p2
            if nurse != curr_nurse
                c2.genotype[patient] = nurse
                construct_solution!(instance, c2, n_nurses)
                route = [instance.patients[p] for p in c2.phenotype[nurse]]
                obj = total_objective(route, instance)
                if obj[1] < best_obj[1] ||
                   (obj[1] == best_obj[1] && obj[2] < best_obj[2]) ||
                   (obj[1] == best_obj[1] && obj[2] == best_obj[2] && obj[3] < best_obj[3])
                    best_nurse = nurse
                    best_obj = obj
                end
            end
        end
        reallocations_c2 = push!(reallocations_c2, best_nurse)
        c2.genotype[patient] = best_nurse
    end

    # println("Patients to reallocate p1: ", nurse_patients_p2)
    # println("Unused nurses: ", unused_nurses_p1)
    # println("Reallocated to: ", reallocations_c1, "\n")
    # println("------------------------------------")

    # println("Patients to reallocate p2: ", nurse_patients_p1)
    # println("Unused nurses: ", unused_nurses_p2)
    # println("Reallocated to: ", reallocations_c2, "\n")
    # println("------------------------------------")


    return c1, c2
end


end # modules