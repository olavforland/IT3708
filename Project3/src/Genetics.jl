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

end

end #module