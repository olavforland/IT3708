# main.jl

include("DataParser.jl")
include("Genetics.jl")
include("Utils.jl")
include("Mutation.jl")
include("Crossover.jl")
include("Selection.jl")
include("TSPHeuristic.jl")
include("VNSHeuristic.jl")
include("GA.jl")

using .DataParser: parse_data, Patient
using .Genetics: Chromosome, compute_fitness!, compute_unfitness!
using .GA: initialize_population, genetic_algorithm
using .Utils: write_chromosome_to_file
using .Mutation: swap_mutation!
using .Crossover: two_point_crossover
using .Selection: tournament_selection, survivor_selection!
# push!(LOAD_PATH, pwd())

instance_nr = 0

# Path to the JSON file with data
readpath = joinpath("data", "train_" * string(instance_nr) * ".json")
writepath = joinpath("solutions", "train_" * string(instance_nr) * ".json")

# Parse the data from the file
instance = parse_data(readpath)

population = genetic_algorithm(instance, 50, 100000, 0.1, instance.n_nurses)
# population = initialize_population(1, instance.n_nurses, instance)

best_individual = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[1]

write_chromosome_to_file(best_individual, writepath)






