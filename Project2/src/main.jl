# main.jl

include("DataParser.jl")
include("Genetics.jl")
include("Utils.jl")
include("GA.jl")
include("Mutation.jl")
include("Crossover.jl")
include("Selection.jl")

using .DataParser: parse_data
using .Genetics: Chromosome, compute_fitness!, compute_unfitness!
using .GA: initialize_population
using .Utils: write_chromosome_to_file
using .Mutation: swap_mutation!
using .Crossover: two_point_crossover
using .Selection: tournament_selection
# push!(LOAD_PATH, pwd())

instance_nr = 0

# Path to the JSON file with data
readpath = joinpath("data", "train_" * string(instance_nr) * ".json")
writepath = joinpath("solutions", "train_" * string(instance_nr) * ".json")

# Parse the data from the file
instance = parse_data(readpath)

test = Chromosome([82, 83, 86, 87, 90, 84, 89, 85, 88, 91])
test.phenotype = [[82, 83, 86, 87, 90, 84, 89, 85, 88, 91]]

compute_fitness!(test, instance)
compute_unfitness!(test, instance)

println(test.fitness)
println(test.time_unfitness)
println(test.strain_unfitness)

swap_mutation!(test)

two_point_crossover(test, test)

population = initialize_population(10, instance.n_nurses, instance)

p1, p2 = tournament_selection(population, 2)



best_individual = population[findmin(getfield.(population, :fitness))[2]]

write_chromosome_to_file(best_individual, writepath)






