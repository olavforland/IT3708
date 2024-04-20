module NSGA2

using ..Genetics: Chromosome


export fast_non_dominated_sort





function fast_non_dominated_sort!(population::Vector{Chromosome})
    """
    Sorts the population into pareto fronts based on edge, deviation and connectivity
    """

    # True if p dominates q
    dominates = (p::Chromosome, q::Chromosome) -> (
        (p.edge <= q.edge && p.deviation <= q.deviation && p.connectivity <= q.connectivity) &&
        (p.edge < q.edge || p.deviation < q.deviation || p.connectivity < q.connectivity)
    )

    pareto_fronts = Vector{Vector{Chromosome}}()
    push!(pareto_fronts, Vector{Chromosome}())

    for p in population
        p.n_dominated = 0
        p.dominated = Set{Chromosome}()

        for q in population
            if dominates(p, q) # If p dominates q
                push!(p.dominated, q) # Add q to the set of solutions dominated by p
            elseif dominates(q, p)
                p.n_dominated += 1 # Increment domination counter of p
            end #if
        end #for
        # p belongs to the first front
        if p.n_dominated == 0
            p.rank = 1
            push!(pareto_fronts[1], p)
        end #if
    end #for

    i = 1
    # While there are solutions in the ith front
    while !isempty(pareto_fronts[i])
        next_front = Vector{Chromosome}()
        for p in pareto_fronts[i]
            # For each solution q dominated by p
            for q in p.dominated
                # Remove domination by p
                q.n_dominated -= 1
                # Check if still dominated
                if q.n_dominated == 0
                    q.rank = i + 1
                    push!(next_front, q)
                end #if
            end #for
        end #for
        i += 1
        if !isempty(next_front)
            push!(pareto_fronts, next_front)
        else
            break
        end

    end #while

    return pareto_fronts
end

function crowding_distance_assignment!(front::Vector{Chromosome})
    """
    Assigns crowding distance to each chromosome in the front
    """

    for p in front
        p.crowding_distance = 0.0
    end #for
    

    # For each objective function
    for objective in [x->x.edge, x->x.connectivity, x->x.deviation]
        front = sort(front, by=objective)

        # Set boundary point distance to inf so they are always selected
        front[1].crowding_distance = Inf
        front[end].crowding_distance = Inf
        
        f_min = objective(front[1])
        f_max = objective(front[end]) 

        # For all other points, add normalized manhatten distance to neighbors
        for i in 2:length(front)-1
            front[i].crowding_distance += (objective(front[i+1]) - objective(front[i-1])) / (f_max - f_min)
        end #for
    end #for
end

function crowded_comparison_operator(i::Chromosome, j::Chromosome)
    """
    Comparison operator for sorting based on crowding distance and rank.
        Individuals in better frontiers are prioritized, with ties broken by crowding distance.
        Returns true if inidividual i is better than individual j
    """
    return (i.rank < j.rank) || ((i.rank == j.rank) && (i.crowding_distance > j.crowding_distance))

end


end # module

