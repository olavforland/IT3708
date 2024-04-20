module MOEA

using Random

using ..Genetics: Chromosome, compute_edge_obj!, compute_connectivity_obj!, compute_deviation_obj!, get_segment_mask
using ..Utils: read_image, mst_to_genotype, min_spanning_tree
using ..Problem: ProblemInstance
using ..NSGA2: fast_non_dominated_sort!, crowding_distance_assignment!, crowded_comparison_operator
using ..Selection: binary_tournament_parent_selection, elitism_survivor_selection
using ..GeneticOperators: uniform_crossover, mutate!

export multi_obj_EA, initialize_population

function multi_obj_EA(population::Vector{Chromosome}, p_mut::Float64, n_generations::Int, n_offspring::Int)
    
    for g in 1:n_generations

        print("Generation ", g, "\n")
        
        # Select parents
        parents = binary_tournament_parent_selection(population, n_offspring)
        
        # Generate offspring
        offspring = Vector{Chromosome}()
        for i in 1:2:length(parents)
            p1, p2 = parents[i:i+1]
            c1, c2 = uniform_crossover(p1, p2)
            mutate!(c1, p_mut), mutate!(c2, p_mut)
            push!(offspring, c1), push!(offspring, c2)
        end #for

        # Select survivors
        population = elitism_survivor_selection(parents, offspring)

    end


end


function initialize_population(n_individuals::Int, instance::ProblemInstance)::Vector{Chromosome}
    """
    Function that initializes a population of n individuals
    """
    println("Initializing population of ", n_individuals, " individuals")
    population = Vector{Chromosome}()

    #create n_inidividuals 
    for _ in 1:n_individuals
        #create mst from pixels and dimensions
        mst = min_spanning_tree(instance.pixels)
        genotype = mst_to_genotype(mst, (instance.height, instance.width))
        chromosome = Chromosome(genotype)
        chromosome.graph = mst
        mask = get_segment_mask(chromosome, instance)
        compute_edge_obj!(chromosome, mask, instance)
        compute_connectivity_obj!(chromosome, mask, instance)
        compute_deviation_obj!(chromosome, mask, instance)
        push!(population, chromosome)

    end #for

    # Run through to calculate rank and crowding distance
    pareto_fronts = fast_non_dominated_sort!(population)
    for front in pareto_fronts
        crowding_distance_assignment!(front)
    end #for

    return population
end

end