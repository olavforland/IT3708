module Similarity

using ..DataParser: ProblemInstance
using ..Genetics: Chromosome

function create_similarity_matrix(population::Vector{Chromosome})
    n = length(population)
    similarity_matrix = zeros(n, n)
    for chromosome_1 in population
        for chromosome_2 in population
            similarity_matrix[chromosome_1.id, chromosome_2.id] = similarity(chromosome_1, chromosome_2)
            similarity_matrix[chromosome_2.id, chromosome_1.id] = similarity(chromosome_2, chromosome_1)
        end
    end
    return similarity_matrix
end

function similarity(chromosome1::Chromosome, chromosome2::Chromosome)
    sim_array = vector{Float64}(undef, length(chromosome1.patient_sets))
    for (i, patient_set) in enumerate(chromosome1.patient_sets)
        sim_array[i] = maximum([len(intersect(patient_set, chromosome2.patient_sets[j])) for j in 1:length(chromosome2.patient_sets)])
    end
    return sum(sim_array)
end

end # module
