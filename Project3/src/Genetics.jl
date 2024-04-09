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
function compute_edge_obj!(chromosome::Chromosome, mask::Vector{Vector{Int}})
    #TODO: 
    #Function for computing edge objective
end

function compute_connectivity_obj!(chromosome::Chromosome, mask::Vector{Vector{Int}})
    #TODO: 
    #Function for computing connectivity objective
end

function compute_deviation_obj!(chrommosome::Chromosome, mask::Vector{Vector{Int}})
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
    segments = Dict{Int,Tuple{Int,Int}}()


    #TODO: Check if it is actually this simple. 
    #Current strat: 
    #1. For each node with children in the forrest, check if it is in a segment.
    #2. If it is in a segment, add the children to the segment.
    #3. If it is not in a segment, create a new segment and add the children to the segment, with the parent node. 
    #4. Repeat until all nodes are in a segment.
    #Could probably be done more efficiently.
    for (r, neighbors) in forrest
        for (r, set) in segments
            if r in set
                push!(set, neighbors)
                break
            end
        end
        segments[segment] = push!(Set{Tuple{Int,Int}}(neighbors), r)
        segment += 1
    end

    #fill mask with segment values
    for (segment, set) in segments
        for (r, c) in set
            mask[r][c] = segment
        end
    end

    return mask
end #get_segment_mask

end #module