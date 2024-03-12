# main.jl

include("DataParser.jl")
include("Genetics.jl")
include("Utils.jl")
include("Mutation.jl")
include("Selection.jl")
include("TSPHeuristic.jl")
include("Objective.jl")
include("LocalSearch.jl")
include("VNSHeuristic.jl")
include("Crossover.jl")
include("LambdaInterchange.jl")
include("LargeNeighborhoodSearch.jl")
include("GA.jl")



using .DataParser: parse_data, Patient
using .Genetics: Chromosome, compute_fitness!, compute_unfitness!
using .GA: initialize_population, genetic_algorithm
using .Utils: write_chromosome_to_file, write_population_to_file
using .Mutation: swap_mutation!
using .Crossover: two_point_crossover
using .Selection: tournament_selection, survivor_selection!
using .LambdaInterchange: lambda_shift_operation
using .LargeNeighborhoodSearch: tsp_all_routes!
# push!(LOAD_PATH, pwd())
using .TSPHeuristic: savelsbergh_heuristic

instance_nr = 2

# Path to the JSON file with data
readpath = joinpath("data", "train_" * string(instance_nr) * ".json")
writepath = joinpath("solutions", "train_" * string(instance_nr) * ".json")

# Parse the data from the file
instance = parse_data(readpath)

# population = initialize_population(30, instance.n_nurses, instance)
population = genetic_algorithm(instance, 30, 50000, 0.3, instance.n_nurses)

# for individual in population
#     tsp_all_routes!(individual, instance)
# end

best_individual = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[1]


write_population_to_file(population, writepath)




#Just testing heuristic.
"""

best_individual.phenotype = [Vector{Int}() for _ in 1:instance.n_nurses]

#for each nurse, generate a solution using the savelsbergh heuristic
for nurse in unique(best_individual.genotype)
    best_individual.phenotype[nurse] = savelsbergh_heuristic(instance, best_individual, nurse)
end

"""