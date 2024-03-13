module LocalSearch

export variable_neighborhood_decent!, local_2_opt!, local_1_shift!, random_1_shift!

using ..Genetics
using ..DataParser


function variable_neighborhood_decent!(route::Vector{Patient}, instance::ProblemInstance, objective::Function)
    prev_route = nothing

    max_iter = 10
    iter = 0

    improved = true

    new_obj = Inf
    
    while improved #&& (iter < max_iter)
        obj = local_1_shift!(route, instance, objective)
        prev_route = map(p -> p.id, route)

        new_obj = local_2_opt!(route, instance, objective)

        iter += 1

        improved = new_obj < obj

    end
    return new_obj

end

function local_2_opt!(route::Vector{Patient}, instance::ProblemInstance, objective::Function)
    best_obj = objective(route, instance)

    for i in 1:length(route)
        for j in i+1:length(route)
            # If j cannot precede i, skip
            if (route[j].id, route[i].id) ∈ instance.inadmissable_presedence
                break
            end

            reverse!(route[i:j])
            
            obj = objective(route, instance)
            if obj < best_obj
                best_obj = obj
            else
                # Change back
                reverse!(route[i:j])
            end
        end
        for j in i-1:-1:1
            # If j cannot precede i, i cannot be inserted after j
            if (route[i].id, route[j].id) ∈ instance.inadmissable_presedence
                break
            end

            reverse!(route[i:j])
            
            obj = objective(route, instance)
            if obj < best_obj
                best_obj = obj
            else
                # Change back
                reverse!(route[i:j])
            end
        end
    end

    return best_obj

end

# ---------------- Helpers ----------------

function local_1_shift!(route::Vector{Patient}, instance::ProblemInstance, objective::Function)
    
    best_obj = objective(route, instance)
    
    # Iterate through all possible 1-shifts
    for i in 1:length(route)
        # Consider the cases when patient i precedes patient j
        for j in i+1:length(route)
            # If j cannot precede i, i cannot be inserted after j
            if (route[j].id, route[i].id) ∈ instance.inadmissable_presedence
                break
            end

            shifted_patient = splice!(route, i) # Remove the element at position i
            insert!(route, j - 1, shifted_patient) # Insert that element at position j
            # Recalculate objective
            obj = objective(route, instance)
            if obj < best_obj
                best_obj = obj
            else
                # Revert the shift operation
                deleteat!(route, j - 1) # Remove the shifted patient from its new position
                insert!(route, i, shifted_patient) # Insert the patient back to its original position
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
            else
                # Revert the shift operation
                deleteat!(route, j) 
                insert!(route, i, shifted_patient) 
            end
        end
    end
    
    return best_obj
end

function random_1_shift!(route::Vector{Patient}, level::Int)

    swaps = Tuple{Int, Int}[]

    # Perform level random 1-shifts 
    for _ in 1:level
        i = rand(1:length(route))
        j = rand(1:length(route))
        if i != j
            route[i], route[j] = route[j], route[i]
            
            # Record the swap
            push!(swaps, (i, j))
        end
    end
    
    return swaps
    
end

end # module