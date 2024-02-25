# main.jl

include("DataParser.jl")
include("Genetics.jl")
include("Modules.jl")
include("Utils.jl")

using .DataParser: parse_data
using .Genetics: Chromosome, compute_fitness!
# push!(LOAD_PATH, pwd())


# Path to the JSON file with data
filepath = joinpath("data", "train_0.json")

# Parse the data from the file
instance = parse_data(filepath)






