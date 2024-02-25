# main.jl
include("DataParser.jl")
using .DataParser

# Path to the JSON file with data
filepath = joinpath("data", "train_0.json")

# Parse the data from the file
instance = parse_data(filepath)

