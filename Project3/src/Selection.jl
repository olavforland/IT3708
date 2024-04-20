

module Selection

using ..Genetics: Chromosome
using ..NSGA2: fast_non_dominated_sort!, crowding_distance_assignment!, crowded_comparison_operator


function binary_tournament_parent_selection(population::Vector{Chromosome}, n_parents::Int=2)::Chromosome
    """
    Function that selects a parent using binary tournament selection
    """

    selected = Vector{Chromosome}()
    for _ in 1:n_parents
        # Randomly select two individuals
        i1 = rand(1:length(population))
        i2 = rand(1:length(population))
        while i2 == i1
            i2 = rand(1:length(population))
        end

        # Select the best individual
        if crowded_comparison_operator(population[i1], population[i2])
            push!(selected, population[i1])
        else
            push!(selected, population[i2])
        end
    end
    return selected
    
end


function elitism_survivor_selection(parents::Vector{Chromosome}, offspring::Vector{Chromosome})

    n = length(parents)

    # Combine parents and offspring
    combined = copy(parents)
    append!(combined, offspring)

    pareto_fronts = fast_non_dominated_sort!(combined)
    next_population = Vector{Chromosome}()
    i = 1

    while length(next_population) + length(pareto_fronts[i]) <= n
        crowding_distance_assignment!(pareto_fronts[i])
        append!(next_population, pareto_fronts[i])
        i += 1
        if i > length(pareto_fronts)
            break
        end
    end

    if length(next_population) < n
        sort!(pareto_fronts[i], by=x->(x.rank, -x.crowding_distance))
        append!(next_population, pareto_fronts[i][1:n-length(next_population)])
    end
    
    return next_population
end

end # module