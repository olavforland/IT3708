module Similarity

using ..DataParser: ProblemInstance
using ..Genetics: Chromosome

export create_similarity_matrix, similarity

function create_similarity_matrix(population::Vector{Chromosome})
    n = length(population)
    similarity_matrix = zeros(n, n)
    for chromosome_1 in population
        for chromosome_2 in population
            avg = (similarity(chromosome_1, chromosome_2) + similarity(chromosome_2, chromosome_1)) / 2
            similarity_matrix[chromosome_1.id, chromosome_2.id], similarity_matrix[chromosome_2.id, chromosome_1.id] = avg, avg
        end
    end
    return similarity_matrix
end

function similarity(chromosome1::Chromosome, chromosome2::Chromosome)
    sim_array = zeros(length(chromosome1.patient_sets))
    for (i, patient_set) in enumerate(chromosome1.patient_sets)
        sim_array[i] = maximum([length(intersect(patient_set, chromosome2.patient_sets[j])) for j in 1:length(chromosome2.patient_sets)])
    end
    return sum(sim_array)
end


function similarity_child_everyone(child::Chromosome, population::Vector{Chromosome})
    child.patient_sets = [Set{Int}(route) for route in child.phenotype]

    child_everyone = zeros(length(population))
    everyone_child = zeros(length(population))
    for (_, chromosome) in enumerate(population)
        child_everyone[chromosome.id] = similarity(child, chromosome)
        everyone_child[chromosome.id] = similarity(chromosome, child)
    end
    return child_everyone, everyone_child
end

function update_similarity_matrix!(similarity_matrix::Matrix{Float64}, individual::Chromosome, population::Vector{Chromosome})
    individual_everyone, everyone_individual = similarity_child_everyone(individual, population)
    similarities = (individual_everyone .+ everyone_individual) / 2
    similarities[individual.id] = length(individual.genotype)
    similarity_matrix[individual.id, :] = similarities
    similarity_matrix[:, individual.id] = similarities
end

end # module

