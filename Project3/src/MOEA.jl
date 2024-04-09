module MOEA

using Random

using ..Genetics: Chromosome, compute_edge_obj!, compute_connectivity_obj!, compute_deviation_obj!
using ..utils: read_image, mst_to_genotype, min_spanning_tree


function multi_obj_GA()

end


function initialize_population(n_individuals::Int, image_path::String)::Vector{Chromosome}
    """
    Function that initializes a population of n individuals
    """
    println("Initializing population of ", n_inidividuals, " individuals")
    population = Vector{Chromosome}()

    # Read image
    pixels = read_image(image_path)
    h = length(pixels)
    w = length(pixels[1])

    #create n_inidividuals 
    for _ in 1:n_individuals
        #create mst from pixels and dimensions
        mst = min_spanning_tree(pixels)
        genotype = mst_to_genotype(mst, (h, w))
        chromosome = Chromosome(genotype)
        compute_edge_obj!(chromosome)
        compute_connectivity_obj!(chromosome)
        compute_deviation_obj!(chromosome)
        chromosome.graph = mst
        chromosome.phenotype = pixels
        push!(population, chromosome)

    end #for
end

end