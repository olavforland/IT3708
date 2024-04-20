# main.jl 


include("Utils.jl")
include("Problem.jl")
include("Genetics.jl")
include("GeneticOperators.jl")
include("EvalInterface.jl")
include("NSGA2.jl")
include("Selection.jl")
include("MOEA.jl")

using .MOEA: multi_obj_EA, initialize_population
using .Genetics: Chromosome
using .Utils: read_image
using .Problem: ProblemInstance
using .EvalInterface: draw_segments
using .NSGA2: fast_non_dominated_sort!, crowding_distance_assignment!, crowded_comparison_operator
using .Selection: binary_tournament_parent_selection, elitism_survivor_selection


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
    multi_obj_EA(population, 0.2, 100, 10)


end #main


main()