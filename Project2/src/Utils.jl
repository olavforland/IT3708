

module Utils

export chromosome_to_dict, write_chromosome_to_file, count_unique_individuals, solution_to_txt

using JSON
using ..Genetics: Chromosome
using ..DataParser: ProblemInstance

function chromosome_to_dict(chromosome::Chromosome)
    return Dict(
        "genotype" => chromosome.genotype,
        "phenotype" => isnothing(chromosome.phenotype) ? nothing : [[int for int in vec] for vec in chromosome.phenotype],
        "fitness" => chromosome.fitness,
        "time_unfitness" => chromosome.time_unfitness,
        "strain_unfitness" => chromosome.strain_unfitness
    )
end

# Serialize and write Chromosome to file
function write_chromosome_to_file(chromosome::Chromosome, filename::String)
    chromo_dict = chromosome_to_dict(chromosome)
    json_str = JSON.json(chromo_dict)
    open(filename, "w") do file
        write(file, json_str)
    end
end

# Serialize and write the entire population to file
function write_population_to_file(population::Vector{Chromosome}, filename::String)
    population_dict = [chromosome_to_dict(chromosome) for chromosome in population]
    json_str = JSON.json(population_dict)
    open(filename, "w") do file
        write(file, json_str)
    end
end


function count_unique_individuals(population::Vector{Chromosome})
    # Convert genotype to hashable string
    genotypes_str = map(c -> join(c.genotype, ","), population)
    
    return length(unique(genotypes_str))
end

function solution_to_txt(best_individual::Chromosome, instance::ProblemInstance, instance_nr::Int)
    # Write the solution to a txt file
    open("solution.txt", "w") do file
        write(file, "Instance: $instance_nr\n")
        write(file, "Fitness: $(best_individual.fitness)\n")
        write(file, "Depot return time: $(instance.depot_return_time) \n\n")
        for (i, route) in enumerate(best_individual.phenotype)
            patients_in_route = [instance.patients[patient] for patient in route]
            if isempty(route)
                write(file, "Nurse $i: Not in use\n")
            else
                total_strain = sum([p.demand for p in patients_in_route])
                travel_time = time_fitness_objective(patients_in_route, instance)
                write(file, "Nurse $i: Covered demand: $total_strain Travel time: $travel_time Route: $(join(route, " "))\n")
            end
        end
    end
end



end # module