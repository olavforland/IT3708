module Genetics

using Statistics
using DataFrames
using ..Utils: euclidean_distance
using ..Problem: ProblemInstance

export Chromosome, compute_edge_obj!, compute_connectivity_obj!, compute_deviation_obj!, get_segment_mask, genotype_to_graph!

mutable struct Chromosome

    #chars corresponding to edges. u, d, l, r, n for up, down, left, right, none
    genotype::Vector{Vector{Char}}


    graph::Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}


    #objectives
    #maxmize value across edges of segments 
    edge::Float64
    #minimize connectivity of segments
    connectivity::Float64
    #minimize deviation of segments mean(segment) - pixel value
    deviation::Float64

    Chromosome(genotype::Vector{Vector{Char}}) = new(genotype, Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}(), 0.0, 0.0, 0.0)
end


#Might be better to combine these to avoid recalculating the same values if they are needed in multiple objectives. 
function compute_edge_obj!(chromosome::Chromosome, mask::Vector{Vector{Int}}, instance::ProblemInstance)
    #TODO: 
    #Function for computing edge objective
    #S.T maximization (originally), but minimized to keep consistent with other objectives. Hence -=.
    h = instance.height
    w = instance.width


    for r in 1:length(h)
        for c in 1:length(w)
            valid_index = (i, j) -> 1 <= i <= length(h) && 1 <= j <= length(w)
            for (di, dj) in [(0, 1), (1, 0), (0, -1), (-1, 0)]
                if valid_index(r + di, c + dj)
                    if mask[r][c] != mask[r+di][c+dj]
                        chromosome.edge -= euclidean_distance(instance.pixels[r][c], instance.pixels[r+di][c+dj])
                    end #if
                end #if
            end #for
        end #for
    end #for
end

function compute_connectivity_obj!(chromosome::Chromosome, mask::Vector{Vector{Int}}, instance::ProblemInstance)
    #TODO: 
    #Function for computing connectivity objective
    #S.T minimization
    for (r, neighbors) in instance.knn
        for (idx, n) in enumerate(neighbors)
            if mask[r[1]][r[2]] != mask[n[1]][n[2]]
                chromosome.connectivity += 1 / idx
            end #if
        end #for
    end #for
end

function compute_deviation_obj!(chromosome::Chromosome, mask::Vector{Vector{Int}}, instance::ProblemInstance)
    #TODO:
    #Function for computing deviation objective
    #S.T minimization
    h = instance.height
    w = instance.width

    num_segments = maximum([maximum(mask[r]) for r in 1:length(h)])


    #mu_k = dict of segment means 0.21932524803251266
    mu_k, segments = get_mu_k(instance.pixels, mask)

    for i in 1:num_segments
        for pixel in segments[i]
            chromosome.deviation += euclidean_distance(pixel, mu_k[i])
        end #for
    end #for
end


function get_segment_mask(chromosome::Chromosome, instance::ProblemInstance)::Vector{Vector{Int}}
    forrest = chromosome.graph
    for (k, v) in forrest
        for tup in v
            if tup in keys(forrest)
                push!(forrest[tup], k)
            else
                forrest[tup] = Set([k])
            end
        end
    end

    h = instance.height
    w = instance.width

    segments = Dict{Int,Set{Tuple{Int,Int}}}()

    v = Set{Tuple{Int,Int}}()

    function bfs!(forrest, root)
        q = [root]
        visited = Set{Tuple{Int,Int}}()
        while !isempty(q)
            node = popfirst!(q)
            if node in visited
                continue
            end
            push!(visited, node)
            if node in keys(forrest)
                for neighbor in forrest[node]
                    push!(q, neighbor)
                end
            end
        end
        v = union!(v, visited)
        return visited
    end

    segmentcounter = 1
    for (k, val) in forrest
        if !(k in v)
            segment = bfs!(forrest, k)
            segments[segmentcounter] = segment
            segmentcounter += 1
        end
    end
    mask = [[0 for _ in 1:w] for _ in 1:h]

    for (k, v) in segments
        for (r, c) in v
            mask[r][c] = k
        end
    end
    return mask
end #get_segment_mask


function genotype_to_graph!(chromosome::Chromosome)
    """
    Function that creates/updates the graph of a chromosome based on the genotype
    """
    genotype = chromosome.genotype
    h = length(chromosome.phenotype)
    w = length(chromosome.phenotype[1])

    graph = Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}()

    valid_index = (i, j) -> 1 <= i <= h && 1 <= j <= w

    for i in 1:h
        for j in 1:w
            if genotype[i][j] == 'n'
                graph[(i, j)] = Set{Tuple{Int,Int}}()
            elseif genotype[i][j] == 'u'
                if valid_index(i - 1, j)
                    if !(i - 1, j) in keys(graph)
                        graph[(i - 1, j)] = Set{Tuple{Int,Int}}()
                    end #if
                    push!(graph[(i - 1, j)], (i, j))
                end #if
            elseif genotype[i][j] == 'd'
                if valid_index(i + 1, j)
                    if !(i + 1, j) in keys(graph)
                        graph[(i + 1, j)] = Set{Tuple{Int,Int}}()
                    end #if
                    push!(graph[(i + 1, j)], (i, j))
                end #if
            elseif genotype[i][j] == 'l'
                if valid_index(i, j - 1)
                    if !(i, j - 1) in keys(graph)
                        graph[(i, j - 1)] = Set{Tuple{Int,Int}}()
                    end #if
                    push!(graph[(i, j - 1)], (i, j))
                end #if
            elseif genotype[i][j] == 'r'
                if valid_index(i, j + 1)
                    if !(i, j + 1) in keys(graph)
                        graph[(i, j + 1)] = Set{Tuple{Int,Int}}()
                    end #if
                    push!(graph[(i, j + 1)], (i, j))
                end #if
            end #if
        end #for
    end #for

    chromosome.graph = graph

end #genotype_to_graph!


function get_mu_k(pixels::Vector{Vector{Tuple{Float64,Float64,Float64}}}, mask::Vector{Vector{Int}})
    """
    Function that calculates the mean of each segment in the mask
    """
    f_pixels = reduce(vcat, pixels)
    f_mask = reduce(vcat, mask)
    zipped = zip(f_pixels, f_mask)

    mu_k = Dict{Int,Tuple{Float64,Float64,Float64}}()
    segments = Dict{Int,Vector{Tuple{Float64,Float64,Float64}}}()

    for (pixel, segment) in zipped
        if segment in keys(segments)
            push!(segments[segment], pixel)
        else
            segments[segment] = [pixel]
        end
    end

    for (k, pixels) in segments
        mu_k[k] = (mean([p[1] for p in pixels]), mean([p[2] for p in pixels]), mean([p[3] for p in pixels]))
    end

    return mu_k, segments

end #

end #module