module MOEA

using Random

using ..Genetics: Chromosome, compute_edge_obj!, compute_connectivity_obj!, compute_deviation_obj!
using ..Utils: read_image, mst_to_genotype, min_spanning_tree
using ..Problem: ProblemInstance

export multi_obj_GA, initialize_population

function multi_obj_GA(population::Vector{Chromosome}, p_cross::Float64, p_mut::Float64, n_generations::Int, n_offspring::Int)
    return population
end


function initialize_population(n_individuals::Int, instance::ProblemInstance)::Vector{Chromosome}
    """
    Function that initializes a population of n individuals
    """
    println("Initializing population of ", n_inidividuals, " individuals")
    population = Vector{Chromosome}()

    #create n_inidividuals 
    for _ in 1:n_individuals
        #create mst from pixels and dimensions
        mst = min_spanning_tree(instance.pixels)
        genotype = mst_to_genotype(mst, (instance.height, instance.width))
        chromosome = Chromosome(genotype)
        chromosome.graph = mst
        compute_edge_obj!(chromosome)
        compute_connectivity_obj!(chromosome)
        compute_deviation_obj!(chromosome)
        push!(population, chromosome)

    end #for
    return population
end

end