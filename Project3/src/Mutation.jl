module Mutation

using ..Genetics: Chromosome

function mutation(chromosome::Chromosome)::Chromosome
    """
    Function that mutates a chromosome

    args: chromosome to be mutated
    returns: mutated chromosome
    """
    idx = rand(1:length(chromosome.genotype))
    direction = rand(['u', 'd', 'l', 'r', 'n'])
    chromosome.genotype[idx] = direction
    return chromosome

end #function

end #module