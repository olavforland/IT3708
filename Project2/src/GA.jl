
module GA

export initialize_population, genetic_algorithm

using ..Crossover: two_point_crossover
using ..VismaCrossover: visma_crossover
using ..DataParser
using ..Genetics
using ..LambdaInterchange: lambda_interchange_operation, lambda_shift_operation
using ..LargeNeighborhoodSearch: tsp_all_routes!
using ..Mutation
using ..VNSHeuristic:
    construct_solution!, improve_single_route, improve_solution!, local_2_opt!
using ..Objective
using ..Selection: survivor_selection!, tournament_selection
using ..TSPHeuristic
using ..Utils: count_unique_individuals
using Printf
using Random




function genetic_algorithm(problem_instance::ProblemInstance, n_individuals::Int, n_generations::Int, mutation_rate::Float64, n_nurses::Int)
    start_time = time()  # Capture start time    

    # Direct initialization
    t0 = time()  # Capture start time
    println("\nInitializing population\n")
    population = initialize_population(n_individuals, n_nurses, problem_instance)
    initialization_time = time() - t0  # Calculate elapsed time
    println("Initialization time: ", @sprintf("%.2f", initialization_time), " s\n")

    construction_time = 0.0
    improvement_time = 0.0
    lambda_shift_time = 0.0
    lambda_interchange_time = 0.0

    for generation in 1:n_generations
        p1, p2 = tournament_selection(population, 2)
        #c1, c2 = two_point_crossover(p1, p2, n_nurses)

        c1, c2 = visma_crossover(p1, p2, n_nurses, problem_instance)

        swap_mutation!(c1, mutation_rate), swap_mutation!(c2, mutation_rate)

        # Perform routing step
        construct_solution!(problem_instance, c1, n_nurses)
        construct_solution!(problem_instance, c2, n_nurses)

        # Evaluate fitness and unfitness
        compute_fitness!(c1, problem_instance), compute_fitness!(c2, problem_instance)
        compute_unfitness!(c1, problem_instance), compute_unfitness!(c2, problem_instance)

        # Get favoured offspring  based on fitness
        child = c1.fitness < c2.fitness ? c1 : c2
        # improvement_time += @elapsed begin
        #     improve_solution!(problem_instance, child)
        # end

        # find average fitness
        avg_fitness = sum(map(c -> c.fitness, population)) / length(population)
        avg_time_unfitness = sum(map(c -> c.time_unfitness, population)) / length(population)
        avg_strain_unfitness = sum(map(c -> c.strain_unfitness, population)) / length(population)


        if generation % 10000 == 0
            best_chromosome = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[1]
            println("Generation: ", generation, " Fitness: ", best_chromosome.fitness, " Time unfitness: ", best_chromosome.time_unfitness, " Strain unfitness: ", best_chromosome.strain_unfitness)
            println("Average fitness: ", avg_fitness, " Average time unfitness: ", avg_time_unfitness, " Average strain unfitness: ", avg_strain_unfitness)

            println("Number of unique individuals: ", count_unique_individuals(population), "/", n_individuals)
            println("Number of nurses used: ", length(unique(best_chromosome.genotype)))
            println("")
        end

        # Check if in population    
        if in(join(child.genotype, ","), map(c -> join(c.genotype, ","), population))
            continue
        end

        if generation % 500000 == 0
            println("\nPerforming large neighborhood search\n")

            # Perform large neighborhood TSP
            for individual in population
                tsp_all_routes!(individual, problem_instance)
            end

            lambda_shift_time += @elapsed begin
                # Perform lambda shift operation
                population = map(c -> lambda_shift_operation(c, problem_instance), population)
            end
            lambda_interchange_time += @elapsed begin
                # Perform lambda interchange operation
                population = map(c -> lambda_interchange_operation(c, problem_instance), population)
            end


            # Finally re-optimize each route with 2-opt procedure
            for chromosome in population
                for i in 1:length(chromosome.phenotype)
                    route = map(p -> problem_instance.patients[p], chromosome.phenotype[i])
                    local_2_opt!(route, problem_instance, total_objective) # Improve
                    chromosome.phenotype[i] = map(p -> p.id, route)
                end
            end


        end

        survivor_selection!(population, child)

    end


    for chromosome in population
        for i in 1:length(chromosome.phenotype)
            route = map(p -> problem_instance.patients[p], chromosome.phenotype[i])
            chromosome.phenotype[i] = map(p -> p.id, route)
        end
    end
    total_time = time() - start_time
    println("\nInitialization time: ", @sprintf("%.2f", initialization_time), " s (", @sprintf("%.2f", initialization_time / total_time * 100), "% of total time)")
    println("Construction time: ", @sprintf("%.2f", construction_time), " s (", @sprintf("%.2f", construction_time / total_time * 100), "% of total time)")
    println("Improvement time: ", @sprintf("%.2f", improvement_time), " s (", @sprintf("%.2f", improvement_time / total_time * 100), "% of total time)")
    println("Lambda shift time: ", @sprintf("%.2f", lambda_shift_time), " s (", @sprintf("%.2f", lambda_shift_time / total_time * 100), "% of total time)")
    println("Lambda interchange time: ", @sprintf("%.2f", lambda_interchange_time), " s (", @sprintf("%.2f", lambda_interchange_time / total_time * 100), "% of total time)")
    println("")


    return population
end


function initialize_population(n_individuals::Int, n_nurses::Int, problem_instance::ProblemInstance)
    population = Vector{Chromosome}()

    n_patients = length(problem_instance.patients)

    patients = sort(problem_instance.patients, by=p -> p.id)
    ranked_patients = sort(problem_instance.patients, by=p -> p.rank)

    # Initial nurse areas
    nurse_areas = create_initial_nurse_areas(n_nurses, n_patients)

    shuffled_patients = shuffle(1:length(problem_instance.patients))
    for i in 1:n_individuals
        delta = shuffled_patients[i]
        chromosome = Chromosome([0 for _ in 1:n_patients], n_nurses)
        chromosome.phenotype = [Vector{Int}() for _ in 1:n_nurses]
        current_nurse = 1

        prev_time_unfitness = 0.0
        prev_strain_unfitness = 0.0

        for j in 1:n_patients

            # Current patient
            p = mod1(j + delta, n_patients)
            patient = ranked_patients[p]

            chromosome.genotype[patient.id] = current_nurse

            construct_solution!(problem_instance, chromosome, n_nurses)

            compute_unfitness!(chromosome, problem_instance)
            compute_fitness!(chromosome, problem_instance)
            # improve_solution!(problem_instance, chromosome)


            if (chromosome.time_unfitness > prev_time_unfitness || chromosome.strain_unfitness > prev_strain_unfitness) #&& (length(chromosome.phenotype[current_nurse]) > n_patients / n_nurses)

                compute_unfitness!(chromosome, problem_instance)
                compute_fitness!(chromosome, problem_instance)

                current_nurse += 1
                current_nurse = mod1(current_nurse, n_nurses)

                prev_time_unfitness = chromosome.time_unfitness
                prev_strain_unfitness = chromosome.strain_unfitness

            end
        end

        # Perform lambda shift operation
        chromosome = lambda_shift_operation(chromosome, problem_instance)
        # Perform lambda interchange operation
        chromosome = lambda_interchange_operation(chromosome, problem_instance)
        # Finally re-optimize each route with 2-opt procedure
        for i in 1:length(chromosome.phenotype)
            route = map(p -> problem_instance.patients[p], chromosome.phenotype[i])
            local_2_opt!(route, problem_instance, total_objective) # Improve
            chromosome.phenotype[i] = map(p -> p.id, route)
        end


        # Add chromosome to population
        compute_fitness!(chromosome, problem_instance)
        compute_unfitness!(chromosome, problem_instance)
        push!(population, chromosome)
    end
    return population
end


# ----------------- Helpers -----------------

# Create initial nurse to ensure that nurses operate roughly in the same area for all individuals in initial population
function create_initial_nurse_areas(N, K)
    # Calculate the number of customers per nurse, and the remainder
    customers_per_nurse = div(K, N)
    remainder = K % N

    # Create the vector with evenly distributed nurse assignments
    assignments = Vector{Int}(undef, K)
    index = 1
    for nurse in 1:N
        # Determine the number of customers for this nurse
        count = customers_per_nurse + (nurse <= remainder ? 1 : 0)
        for _ in 1:count
            assignments[index] = nurse
            index += 1
        end
    end

    return assignments
end

end # module

