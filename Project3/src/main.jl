# main.jl 


include("Utils.jl")
include("Problem.jl")
include("Genetics.jl")
include("MOEA.jl")
include("GeneticOperators.jl")
include("EvalInterface.jl")


using .MOEA: multi_obj_GA, initialize_population
using .Genetics: Chromosome
using .Utils: read_image
using .Problem: ProblemInstance
using .EvalInterface: draw_segments

function main()
    print("Main function\n")
    image_number = 86016
    image_path = "images/" * string(image_number) * "/Test image.jpg"

    probleminstance = ProblemInstance(image_path, 8)
    #args: n_individuals, image_path
    population = initialize_population(10, probleminstance)

    for i in 1:10
        println(population[i].edge, " ", population[i].connectivity, " ", population[i].deviation)
    end #for

    draw_segments(population[1], probleminstance)

    # args: population, p_cross, p_mut, n_generations, n_offspring
    multi_obj_GA(population, 0.8, 0.1, 100, 10)



end #main


main()