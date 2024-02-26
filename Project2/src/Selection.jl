module Selection

export tournament_selection

using ..Genetics: Chromosome


function tournament_selection(population::Vector{Chromosome}, n::Int)
    selected = Vector{Chromosome}()
    for _ in 1:n
        # Randomly select two individuals
        i1 = rand(1:length(population))
        i2 = rand(1:length(population))
        while i2 == i1
            i2 = rand(1:length(population))
        end

        # Select the best individual
        if population[i1].fitness > population[i2].fitness
            push!(selected, population[i1])
        else
            push!(selected, population[i2])
        end
    end
    return selected
end

end # module