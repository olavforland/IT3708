module Genetics

export Chromosome, compute_edge_obj!, compute_connectivity_obj!, compute_deviation_obj!

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

function compute_connectivity_obj!(chromosome::Chromosome)
    #TODO: 
    #Function for computing connectivity objective
end

function compute_deviation_obj!(chrommosome::Chromosome)
    #TODO:
    #Function for computing deviation objective
end


function get_segment_mask(chromosome::Chromosome)::Vector{Vector{Int}}
    """
    Function that returns a mask of the segments in the chromosome
    Segments are represented by integers. 
    [[1,1,1,2,2,2],
     [1,1,1,2,2,2],
     [1,1,1,2,2,2],
     [3,3,3,4,4,4],
     [3,3,3,4,4,4],
     [3,3,3,4,4,4]]
    """

    forrest = chromosome.graph
    h = length(chromosome.phenotype)
    w = length(chromosome.phenotype[1])

    #initialize mask
    mask = Vector{Vector{Int}}(fill(0, h, w))

    #initialize segment counter
    segment = 1

    #initialize dictionary to keep track of segments
    segments = Dict{Tuple{Int,Int},Int}()


    #TODO: Check if it is actually this simple. 
    for (r, neighbors) in forrest
        for v in neighbors
            if haskey(segments, r)
                segments[v] = segments[r]
            else
                segments[r] = segment
                segments[v] = segment
                segment += 1
            end #if
        end #for
    end #for

end #get_segment_mask

end #module