module Crossover

export two_point_crossover

using ..Genetics: Chromosome

function two_point_crossover(p1::Chromosome, p2::Chromosome)

    n_patients = length(p1.genotype)

    # Randomly select two points
    point1 = rand(1:n_patients)
    point2 = rand(1:n_patients)

    # Ensure point1 < point2
    if point1 > point2
        point1, point2 = point2, point1
    end

        
    # Create children
    c1 = Chromosome([0 for _ in 1:n_patients])
    c2 = Chromosome([0 for _ in 1:n_patients])


    # Copy the genes from the parents
    for i in 1:n_patients
        if i < point1 || i > point2
            c1.genotype[i] = p1.genotype[i]
            c2.genotype[i] = p2.genotype[i]
        else
            c1.genotype[i] = p2.genotype[i]
            c2.genotype[i] = p1.genotype[i]
        end
    end


    return c1, c2

end


end # module