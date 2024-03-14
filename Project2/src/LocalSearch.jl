module LocalSearch

export variable_neighborhood_decent!, local_2_opt!, local_3_opt!, local_1_shift!, random_1_shift!

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
    improved = true
    

    while improved
        improved = false
        for i in 1:length(route)
            for j in i+1:length(route)
                # If j cannot precede i, skip
                if (route[j].id, route[i].id) ∈ instance.inadmissable_presedence
                    break
                end

                reverse!(route, i, j)
                # route[i:j] = reverse(route[i:j])
                
                obj = objective(route, instance)
                if obj < best_obj
                    best_obj = obj
                    improved = true
                    break
                else
                    # Change back
                    reverse!(route, i, j)
                    # route[i:j] = reverse(route[i:j])
                end
            end
            if improved
                break
            end
        end

    end
    return best_obj
end

function local_3_opt!(route::Vector{Patient}, instance::ProblemInstance, objective::Function)
    best_obj = objective(route, instance)
    n = length(route)
    improved = true
    while improved
        improved = false
        for c1 in 1:n
            i = c1
            for c2 in 2:n-2
                j = mod1(i + c2, n)
                
                for c3 in c2+1:n

                    k = mod1(i + c3, n)

                    # If j cannot precede i, skip
                    if (route[mod1(j+1, n)].id, route[k].id) ∉ instance.inadmissable_presedence &&
                        (route[mod1(i+1, n)].id, route[j].id) ∉ instance.inadmissable_presedence

                        obj = reverse_double_segment!(route, mod1(j+1, n), k, mod1(i+1, n), j, best_obj, instance, objective)
                        if obj < best_obj
                            # println("Improved objective from $best_obj to $obj")
                            best_obj = obj
                            improved = true
                            break
                        end
                    end

                    
                    if (route[j].id, route[mod1(k+1, n)].id) ∉ instance.inadmissable_presedence &&
                        (route[j].id, route[i].id) ∉ instance.inadmissable_presedence &&
                        (route[j].id, route[mod1(i+1, n)].id) ∉ instance.inadmissable_presedence

                        obj = reverse_double_segment!(route, mod1(k+1, n), i, mod1(i+1, n), j, best_obj, instance, objective)
                    
                        if obj < best_obj
                            # println("Improved objective from $best_obj to $obj")
                            best_obj = obj
                            improved = true
                            break
                        end
                    end

                    if (route[mod1(j+1, n)].id, route[mod1(k+1, n)].id) ∉ instance.inadmissable_presedence &&
                        (route[mod1(j+1, n)].id, route[i].id) ∉ instance.inadmissable_presedence &&
                        (route[k].id, route[mod1(k+1, n)].id) ∉ instance.inadmissable_presedence &&
                        (route[k].id, route[i].id) ∉ instance.inadmissable_presedence 

                        if (route[i].id, route[mod1(k+1, n)].id) ∉ instance.inadmissable_presedence
                            
                            obj = reverse_double_segment!(route, mod1(k+1, n), i, mod1(j+1, n), k, best_obj, instance, objective)
                            if obj < best_obj  
                                # println("Improved objective from $best_obj to $obj") 
                                best_obj = obj
                                improved = true
                                break
                            end
                        end
    
                        reverse!(route, mod1(k+1, n), i)
                        reverse!(route, mod1(i+1, n), j)
                        reverse!(route, mod1(j+1, n), k)
    
                        obj = objective(route, instance)
                        if obj < best_obj
                            # println("Improved objective from $best_obj to $obj")
                            best_obj = obj
                            improved = true
                            break
                        else
                            # Change back
                            reverse!(route, mod1(j+1, n), k)
                            reverse!(route, mod1(i+1, n), j)
                            reverse!(route, mod1(k+1, n), i)
                        end

                    end 
                end
                if improved
                    break
                end
            end
            if improved
                break
            end
        end
    end

end


# ---------------- Helpers ----------------

function reverse_single_segment!(route::Vector{Patient}, i::Int, j::Int, best_obj::Tuple{Float64, Float64, Float64}, instance::ProblemInstance, objective::Function)

    reverse!(route, i, j)

    obj = objective(route, instance)

    if obj < best_obj
        return obj
    end
    # Change back
    reverse!(route, i, j)
    return best_obj
end

function reverse_double_segment!(route::Vector{Patient}, i::Int, j::Int, k::Int, l::Int, best_obj::Tuple{Float64, Float64, Float64}, instance::ProblemInstance, objective::Function)

    reverse!(route, i, j)
    reverse!(route, k, l)

    obj = objective(route, instance)

    if obj < best_obj
        return obj
    end
    # Change back
    reverse!(route, k, l)
    reverse!(route, i, j)
    
    return best_obj
end


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