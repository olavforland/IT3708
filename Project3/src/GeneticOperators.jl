module GeneticOperators

using ..Genetics: Chromosome


export crossover, mutation

function single_point_crossover(p1::Chromosome, p2::Chromosome)::Tuple{Chromosome,Chromosome}
    """
    Function that performs single point crossover between two chromosomes

    args: two chromosomes to be crossed over
    returns: two new chromosomes
    """
    idx = rand(1:length(p1.genotype))
    c1 = Chromosome(p1.genotype[1:idx] + p2.genotype[idx+1:end])
    c2 = Chromosome(p2.genotype[1:idx] + p1.genotype[idx+1:end])
    return (c1, c2)

end


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