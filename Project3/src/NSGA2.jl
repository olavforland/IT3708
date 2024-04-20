module NSGA2

using ..Genetics: Chromosome


export fast_non_dominated_sort





function fast_non_dominated_sort(population::Vector{Chromosome})
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
        push!(pareto_fronts, next_front)
    end #while

    return pareto_fronts
end



end # module