module GeneticOperators

using ..Genetics: Chromosome


export uniform_crossover, mutation


function parent_selection(population::Vector{Chromosome})
    """
    """
end




function uniform_crossover(p1::Chromosome, p2::Chromosome)::Tuple{Chromosome,Chromosome}
    #TODO:
    #Implement uniform crossover
    c1 = deepcopy(p1)
    c2 = deepcopy(p2)

    mask = [rand(Bool, length(p1.genotype[1])) for _ in 1:length(p1.genotype)]

    for i in 1:length(p1.genotype)
        for j in 1:length(p1.genotype[1])
            if mask[i][j]
                c1.genotype[i][j] = p2.genotype[i][j]
                #Update graph dict
                update_graph_dict!(c1, (i, j), p1.genotype[i][j], p2.genotype[i][j])
                c2.genotype[i][j] = p1.genotype[i][j]
                #Update graph dict
                update_graph_dict!(c2, (i, j), p2.genotype[i][j], p1.genotype[i][j])
            end #if
        end #for
    end #for

    return (c1, c2)

end


function mutation(chromosome::Chromosome)::Chromosome
    """
    Function that mutates a chromosome

    args: chromosome to be mutated
    returns: mutated chromosome
    """
    i = rand(1:length(chromosome.genotype))
    j = rand(1:length(chromosome.genotype[1]))
    direction = rand(['u', 'd', 'l', 'r', 'n'])
    old = chromosome.genotype[i][j]
    new = direction
    chromosome.genotype[i][j] = direction
    update_graph_dict!(chromosome, (i, j), old, new)
    return chromosome

end #function

function update_graph_dict!(chromosome::Chromosome, changed_idx::Tuple{Int,Int}, old::Char, new::Char)
    """
    Function that updates the segment dictionary of a chromosome after a mutation

    args: chromosome, index of the changed gene, old direction, new direction
    returns: updated segment dictionary
    """
    char_to_dir = Dict('u' => (-1, 0), 'd' => (1, 0), 'l' => (0, -1), 'r' => (0, 1), 'n' => (0, 0))

    if old == new
        return
    end #if

    if (old == 'n' && new != "n")
        dir = char_to_dir[new]
        new_pointed_to = (changed_idx[1] + dir[1], changed_idx[2] + dir[2])
        if !haskey(chromosome.graph, new_pointed_to)
            chromosome.graph[new_pointed_to] = Set{Tuple{Int,Int}}()
        end #if
        push!(chromosome.graph[new_pointed_to], changed_idx)
    end #if

    if (old != 'n' && new == "n")
        dir = char_to_dir[old]
        old_pointed_to = (changed_idx[1] + dir[1], changed_idx[2] + dir[2])
        delete!(chromosome.graph[old_pointed_to], changed_idx)
        if isempty(chromosome.graph[old_pointed_to])
            delete!(chromosome.graph, old_pointed_to)
        end #if
        if !haskey(chromosome.graph, changed_idx)
            chromosome.graph[changed_idx] = Set{Tuple{Int,Int}}()
        end #if
    end #if

    if (old != 'n' && new != 'n')
        new_dir = char_to_dir[new]
        old_dir = char_to_dir[old]
        new_pointed_to = (changed_idx[1] + new_dir[1], changed_idx[2] + new_dir[2])
        if !haskey(chromosome.graph, new_pointed_to)
            chromosome.graph[new_pointed_to] = Set{Tuple{Int,Int}}()
        end #if
        push!(chromosome.graph[new_pointed_to], changed_idx)
        old_pointed_to = (changed_idx[1] + old_dir[1], changed_idx[2] + old_dir[2])
        delete!(chromosome.graph[old_pointed_to], changed_idx)
        if isempty(chromosome.graph[old_pointed_to])
            delete!(chromosome.graph, old_pointed_to)
        end #if
    end #if

end #function


end #module