# main.jl

include("DataParser.jl")
include("Genetics.jl")
include("Utils.jl")
include("Mutation.jl")
include("Crossover.jl")
include("Selection.jl")
include("GA.jl")

using .DataParser: parse_data, Patient
using .Genetics: Chromosome, compute_fitness!, compute_unfitness!
using .GA: initialize_population, genetic_algorithm
using .Utils: write_chromosome_to_file
using .Mutation: swap_mutation!
using .Crossover: two_point_crossover
using .Selection: tournament_selection, partition_population, survivor_selection!
# push!(LOAD_PATH, pwd())

instance_nr = 0

# Path to the JSON file with data
readpath = joinpath("data", "train_" * string(instance_nr) * ".json")
writepath = joinpath("solutions", "train_" * string(instance_nr) * ".json")

# Parse the data from the file
instance = parse_data(readpath)


population = genetic_algorithm(instance, 50, 100000000, 0.1, instance.n_nurses)

best_individual = population[findmin(getfield.(population, :fitness))[2]]

write_chromosome_to_file(best_individual, writepath)






