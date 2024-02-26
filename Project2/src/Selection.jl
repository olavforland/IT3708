module Selection

export tournament_selection, partition_population

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

function survivor_selection!(population::Vector{Chromosome}, child::Chromosome)
    subsets = partition_population(population, child)

    for subset in subsets
        if !isempty(subset)
            # Sort the subset based on your criteria
            sort!(subset, by=x -> (x.time_unfitness, x.strain_unfitness, x.fitness), rev=true)
            # The worst chromosome is now the first in the sorted subset
            worst_chromosome = subset[1]
            println("Worst chromosome: ", worst_chromosome)
            println("")
            println("Child: ", child)


            # Find the index of this chromosome in the population
            worst_index = findfirst(x -> x == worst_chromosome, population)
            
            if worst_index !== nothing
                population[worst_index] = child
                break
            end
        end
    end
end

# ------------------ Helpers ------------------ #
function partition_population(population::Vector{Chromosome}, ref::Chromosome)
    # Initialize subsets
    subsets = [Vector{Chromosome}() for _ in 1:7]

    for individual in population
        # Determine the subset based on conditions
        if (individual.fitness >= ref.fitness) && (individual.strain_unfitness >= ref.strain_unfitness) && (individual.time_unfitness >= ref.time_unfitness)
            push!(subsets[1], individual)
        elseif (individual.fitness < ref.fitness) && (individual.strain_unfitness >= ref.strain_unfitness) && (individual.time_unfitness >= ref.time_unfitness)
            push!(subsets[2], individual)
        elseif (individual.fitness >= ref.fitness) && (individual.strain_unfitness < ref.strain_unfitness) && (individual.time_unfitness >= ref.time_unfitness)
            push!(subsets[3], individual)
        elseif (individual.fitness < ref.fitness) && (individual.strain_unfitness < ref.strain_unfitness) && (individual.time_unfitness >= ref.time_unfitness)
            push!(subsets[4], individual)
        elseif (individual.fitness >= ref.fitness) && (individual.strain_unfitness >= ref.strain_unfitness) && (individual.time_unfitness < ref.time_unfitness)
            push!(subsets[5], individual)
        elseif (individual.fitness < ref.fitness) && (individual.strain_unfitness >= ref.strain_unfitness) && (individual.time_unfitness < ref.time_unfitness)
            push!(subsets[6], individual)
        elseif (individual.fitness >= ref.fitness) && (individual.strain_unfitness < ref.strain_unfitness) && (individual.time_unfitness < ref.time_unfitness)
            push!(subsets[7], individual)
        # else
        #     push!(subsets[8], individual)
        end
    end
    return subsets
end

end # module