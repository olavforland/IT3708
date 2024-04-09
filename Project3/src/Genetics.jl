module Genetics

mutable struct Chromosome

    #chars corresponding to edges. u, d, l, r, n for up, down, left, right, none
    genotype::Vector{Char}

    #vector of vectors of tuples of RGB values
    phenotype::Vector{Vector{Tuple{float64,float64,float64}}}

    graph::Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}


    #objectives
    #maxmize value across edges of segments 
    edge::Float64
    #minimize connectivity of segments
    connectivity::Float64
    #minimize deviation of segments mean(segment) - pixel value
    deviation::Float64

    Chromosome(genotype::Vector{Char}) = new(genotype, nothing, nothing, nothing, nothing, nothing)
end


#Might be better to combine these to avoid recalculating the same values if they are needed in multiple objectives. 
function compute_edge_obj!(chromosome::Chromosome)
    #TODO: 
    #Function for computing edge objective
end

function compute_connectivity_obj(chromosome::Chromosome)
    #TODO: 
    #Function for computing connectivity objective
end

function compute_deviation_obj(chrommosome::Chromosome)
    #TODO:
    #Function for computing deviation objective
end

end #module