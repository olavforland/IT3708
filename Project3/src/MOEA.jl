module MOEA

using Random
using Base.Threads

using ..Genetics: Chromosome, compute_edge_obj!, compute_connectivity_obj!, compute_deviation_obj!, get_segment_mask
using ..Utils: read_image, mst_to_genotype, min_spanning_tree
using ..Problem: ProblemInstance
using ..NSGA2: fast_non_dominated_sort!, crowding_distance_assignment!, crowded_comparison_operator
using ..Selection: binary_tournament_parent_selection, elitism_survivor_selection
using ..GeneticOperators: uniform_crossover, mutate!


export multi_obj_EA, initialize_population

function multi_obj_EA(population::Vector{Chromosome}, p_mut::Float64, n_generations::Int, n_offspring::Int, instance::ProblemInstance)
    
    for g in 1:n_generations

        print("-----------Generation ", g, "-----------\n")
        
        # Select parents
        parents = binary_tournament_parent_selection(population, n_offspring)
        
        # Generate offspring
        offspring = Vector{Chromosome}()
        @threads for i in 1:2:length(parents)
            p1, p2 = parents[i:i+1]
            c1, c2 = uniform_crossover(p1, p2)
            mutate!(c1, p_mut), mutate!(c2, p_mut)

            mask_c1 = get_segment_mask(c1, instance)
            mask_c2 = get_segment_mask(c2, instance)

            compute_edge_obj!(c1, mask_c1, instance), compute_edge_obj!(c2, mask_c2, instance)
            compute_connectivity_obj!(c1, mask_c1, instance), compute_connectivity_obj!(c2, mask_c2, instance)
            compute_deviation_obj!(c1, mask_c1, instance), compute_deviation_obj!(c2, mask_c2, instance)

            push!(offspring, c1), push!(offspring, c2)

            # println("Offspring ", i, " edge obj: ", c1.edge, " connectivity obj: ", c1.connectivity, " deviation obj: ", c1.deviation)
            # println("Offspring ", i+1, " edge obj: ", c2.edge, " connectivity obj: ", c2.connectivity, " deviation obj: ", c2.deviation)

            s = 0
            for i in 1:length(c1.genotype)
                
                s += length(findall(x->x=='n', c1.genotype[i]))

            end
            # println("Number of none in offspring ", i, ": ", s)
        end #for

        # Select survivors
        population = elitism_survivor_selection(parents, offspring)

        sorted_pop = sort(population, by=x->(x.edge, x.connectivity, x.deviation))
        println("Best individual in generation ", g, " has edge obj: ", sorted_pop[1].edge, " connectivity obj: ", sorted_pop[1].connectivity, " deviation obj: ", sorted_pop[1].deviation) 
        println("Worst individual in generation ", g, " has edge obj: ", sorted_pop[end].edge, " connectivity obj: ", sorted_pop[end].connectivity, " deviation obj: ", sorted_pop[end].deviation) 

    end

    return population


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