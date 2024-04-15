# main.jl 


include("Utils.jl")
include("Problem.jl")
include("Genetics.jl")
include("MOEA.jl")
include("GeneticOperators.jl")


using .MOEA
using .Genetics
using .Utils
using .Problem: ProblemInstance

function main()
    image_number = 86016
    image_path = "images/" * string(image_number) * "/Test image.jpg"

    probleminstance = ProblemInstance(image_path, 8)

    #args: n_individuals, image_path
    population = initialize_population(10, probleminstance)

    # args: population, p_cross, p_mut, n_generations, n_offspring
    multi_obj_GA(population, 0.8, 0.1, 100, 10)



end #main
