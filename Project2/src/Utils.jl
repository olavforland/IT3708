

module Utils

export chromosome_to_dict, write_chromosome_to_file, count_unique_individuals

using JSON
using ..Genetics: Chromosome


function chromosome_to_dict(chromosome::Chromosome)
    return Dict(
        "genotype" => chromosome.genotype,
        "phenotype" => isnothing(chromosome.phenotype) ? nothing : [[int for int in vec] for vec in chromosome.phenotype],
        "fitness" => chromosome.fitness,
        "time_unfitness" => chromosome.time_unfitness,
        "strain_unfitness" => chromosome.strain_unfitness
    )
end

# Serialize and write Chromosome to file
function write_chromosome_to_file(chromosome::Chromosome, filename::String)
    chromo_dict = chromosome_to_dict(chromosome)
    json_str = JSON.json(chromo_dict)
    open(filename, "w") do file
        write(file, json_str)
    end
end


function count_unique_individuals(population::Vector{Chromosome})
    # Convert genotype to hashable string
    genotypes_str = map(c -> join(c.genotype, ","), population)
    
    return length(unique(genotypes_str))
end

end # module