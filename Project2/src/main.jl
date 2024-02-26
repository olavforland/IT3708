# main.jl

include("DataParser.jl")
include("Genetics.jl")
include("Utils.jl")

using .DataParser: parse_data
using .Genetics: Chromosome, compute_fitness!, compute_unfitness!
# push!(LOAD_PATH, pwd())


# Path to the JSON file with data
filepath = joinpath("data", "train_0.json")

# Parse the data from the file
instance = parse_data(filepath)

test = Chromosome([82, 83, 86, 87, 90, 84, 89, 85, 88, 91])
test.phenotype = [[82, 83, 86, 87, 90, 84, 89, 85, 88, 91]]

compute_fitness!(test, instance)
compute_unfitness!(test, instance)

println(test.fitness)
println(test.time_unfitness)
println(test.strain_unfitness)





