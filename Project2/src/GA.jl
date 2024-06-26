
module GA

export initialize_population, genetic_algorithm

using Random
using Printf
using Base.Threads
using Distributed

using ..DataParser
using ..Genetics
using ..Similarity: create_similarity_matrix, similarity_child_everyone
using ..Crossover
using ..Mutation
using ..Selection: tournament_selection, similarity_selection, survivor_selection!
using ..TSPHeuristic
using ..VNSHeuristic: construct_solution!, construct_single_route, improve_solution!, local_2_opt!, local_3_opt!, improve_single_route
using ..Utils: count_unique_individuals
using ..LambdaInterchange: lambda_shift_operation, lambda_interchange_operation
using ..Objective
using ..LargeNeighborhoodSearch: tsp_all_routes!



function genetic_algorithm(initial_population::Vector{Chromosome}, problem_instance::ProblemInstance, n_individuals::Int, n_generations::Int, mutation_rate::Float64, crossover_op::Function, n_nurses::Int, lns_frequency::Int=10000, verbose::Int=1, print_frequency::Int=1000)
    start_time = time()  # Capture start time    

    population = deepcopy(initial_population)

    improvement_time = 0.0
    lambda_shift_time = 0.0
    lambda_interchange_time = 0.0

    # similarity_matrix = create_similarity_matrix(population)

    for generation in 1:n_generations

        n_offspring = 7 * n_individuals

        offspring = Vector{Chromosome}()


        # Generate n_offspring children
        @threads for i in 1:n_offspring


            p1, p2 = tournament_selection(population, 2)
            # Use visma crossover 15% of the time, otherwise randomly select uniform, 1-point, two-point or three-point crossover
            if rand() < 0.15
                c1, c2 = visma_crossover(p1, p2, n_nurses, problem_instance)
            else
                c1, c2 = crossover_op(rand([0, 1, 2, 3]))(p1, p2, problem_instance)
            end

            # Mutate
            swap_mutation!(c1, mutation_rate), swap_mutation!(c2, mutation_rate)
            insert_mutation!(c1, mutation_rate), insert_mutation!(c2, mutation_rate)

            # Perform routing step
            construct_solution!(problem_instance, c1, n_nurses)
            construct_solution!(problem_instance, c2, n_nurses)

            # Improve routes
            for j in 1:length(c1.phenotype)
                route = map(p -> problem_instance.patients[p], c1.phenotype[j])
                local_2_opt!(route, problem_instance, total_objective) # Improve
                # local_3_opt!(route, problem_instance, total_objective) # Improve
                c1.phenotype[j] = map(p -> p.id, route)
            end

            for j in 1:length(c2.phenotype)
                route = map(p -> problem_instance.patients[p], c2.phenotype[j])
                local_2_opt!(route, problem_instance, total_objective) # Improve
                # local_3_opt!(route, problem_instance, total_objective) # Improve
                c2.phenotype[j] = map(p -> p.id, route)
            end


            # Evaluate fitness and unfitness
            compute_fitness!(c1, problem_instance), compute_fitness!(c2, problem_instance)
            compute_unfitness!(c1, problem_instance), compute_unfitness!(c2, problem_instance)

            # Get favoured offspring  based on fitness
            child = c1.fitness < c2.fitness ? c1 : c2

            # Push to offspring
            push!(offspring, child)
        end



        if generation % print_frequency == 0 && verbose > 0
            best_chromosome = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[1]
            top5_chromosome = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[5]
            top10_chromosome = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[10]
            worst_chromosome = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[end]
            println("-------------- Generation: ", generation, "-------------- ")
            println("Top 1 - Fitness: ", @sprintf("%.2f", best_chromosome.fitness), " Time unfitness: ", @sprintf("%.2f", best_chromosome.time_unfitness), " Strain unfitness: ", @sprintf("%.2f", best_chromosome.strain_unfitness))
            println("Top 5 - Fitness: ", @sprintf("%.2f", top5_chromosome.fitness), " Time unfitness: ", @sprintf("%.2f", top5_chromosome.time_unfitness), " Strain unfitness: ", @sprintf("%.2f", top5_chromosome.strain_unfitness))
            println("Top 10 - Fitness: ", @sprintf("%.2f", top10_chromosome.fitness), " Time unfitness: ", @sprintf("%.2f", top10_chromosome.time_unfitness), " Strain unfitness: ", @sprintf("%.2f", top10_chromosome.strain_unfitness))
            println("Top 30 - Fitness: ", @sprintf("%.2f", worst_chromosome.fitness), " Time unfitness: ", @sprintf("%.2f", worst_chromosome.time_unfitness), " Strain unfitness: ", @sprintf("%.2f", worst_chromosome.strain_unfitness))
            println("Number of unique individuals: ", count_unique_individuals(population), "/", n_individuals)
            println("Number of nurses used: ", length(unique(best_chromosome.genotype)))
            println("---------------------------------------------\n")
        end



        if generation % lns_frequency == 0
            if verbose > 0
                println("\nPerforming large neighborhood search\n")
            end


            @threads for i in eachindex(population)
                # Perform large neighborhood TSP
                population[i] = tsp_all_routes!(population[i], problem_instance)

                lambda_shift_time += @elapsed begin
                    population[i] = lambda_shift_operation(population[i], problem_instance)
                end
                lambda_interchange_time += @elapsed begin
                    population[i] = lambda_interchange_operation(population[i], problem_instance)
                end

                # Finally re-optimize each route with 2-opt procedure
                for j in 1:length(population[i].phenotype)
                    route = map(p -> problem_instance.patients[p], population[i].phenotype[j])
                    local_2_opt!(route, problem_instance, total_objective) # Improve
                    population[i].phenotype[j] = map(p -> p.id, route)
                end
            end

            # similarity_matrix = create_similarity_matrix(population)

        end
        # # Check if in population    
        # if join(child.genotype, ",") ∉ map(c -> join(c.genotype, ","), population)
        #     survivor_selection!(population, child, similarity_matrix)
        # end


        # Make similar offspring compete
        similarity_matrix = create_similarity_matrix(offspring)
        while length(offspring) > n_individuals * 3

            # Randomly select individual
            ind = rand(1:length(offspring))
            # Select most similar individual
            min_ind = argmax(similarity_matrix[1:end.!=ind, ind])

            # Remove both children from offspring
            c1 = splice!(offspring, ind)
            c2 = splice!(offspring, ind < min_ind ? min_ind - 1 : min_ind)

            # Pick survivor
            child = c1.fitness < c2.fitness ? c1 : c2
            dead_ind = c1.fitness < c2.fitness ? min_ind : ind

            push!(offspring, child)

            # Remove row and column from similarity matrix for dead child
            similarity_matrix = similarity_matrix[1:end.!=dead_ind, 1:end.!=dead_ind]

        end
        # Take the n_individuals best individuals
        # population = sort(offspring, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[1:n_individuals]

        # Select survivors
        for child in offspring
            if join(child.genotype, ",") ∉ map(c -> join(c.genotype, ","), population)
                survivor_selection!(population, child, similarity_matrix)
            end
        end

        # population = offspring


    end


    for chromosome in population
        for i in 1:length(chromosome.phenotype)
            route = map(p -> problem_instance.patients[p], chromosome.phenotype[i])
            local_2_opt!(route, problem_instance, total_objective) # Improve
            chromosome.phenotype[i] = map(p -> p.id, route)
        end
    end



    return population
end


function initialize_population(n_individuals::Int, n_nurses::Int, problem_instance::ProblemInstance, verbose::Int=1)

    if verbose > 0
        println("\nInitializing population with " * string(nthreads()) * " threads\n")
    end

    population = Vector{Chromosome}()

    n_patients = length(problem_instance.patients)

    patients = sort(problem_instance.patients, by=p -> p.id)
    ranked_patients = sort(problem_instance.patients, by=p -> p.rank)
    # Initial nurse areas
    nurse_areas = create_initial_nurse_areas(n_nurses, n_patients)

    shuffled_patients = shuffle(1:length(patients))

    @threads for i in 1:n_individuals
        delta = shuffled_patients[i]
        current_nurse = nurse_areas[delta]

        chromosome = Chromosome([0 for _ in 1:n_patients], n_nurses)
        chromosome.phenotype = [Vector{Int}() for _ in 1:n_nurses]

        for j in 1:n_patients

            # Current patient
            p = mod1(j + delta, n_patients)
            patient = ranked_patients[p]

            chromosome.genotype[patient.id] = current_nurse

            # construct_solution!(problem_instance, chromosome, n_nurses)

            compute_unfitness!(chromosome, problem_instance)
            compute_fitness!(chromosome, problem_instance)

            # Improve route
            for k in 1:length(chromosome.phenotype)
                route = map(p -> patients[p], chromosome.phenotype[k])
                local_2_opt!(route, problem_instance, total_objective) # Improve
                chromosome.phenotype[k] = map(p -> p.id, route)
            end


            if mod1(j, 5) == 20 #(chromosome.time_unfitness > 0.0 || chromosome.strain_unfitness > 0.0) #&& (length(chromosome.phenotype[current_nurse]) > n_patients / n_nurses)

                construct_solution!(problem_instance, chromosome, n_nurses)

                current_nurse += 1
                current_nurse = mod1(current_nurse, n_nurses)



                # Recompute (un)fitness
                compute_unfitness!(chromosome, problem_instance)
                compute_fitness!(chromosome, problem_instance)
            end
        end


        # Perform lambda shift operation
        chromosome = lambda_shift_operation(chromosome, problem_instance)
        # Perform lambda interchange operation
        chromosome = lambda_interchange_operation(chromosome, problem_instance)
        # Finally re-optimize each route with 2-opt procedure
        for i in 1:length(chromosome.phenotype)
            route = map(p -> patients[p], chromosome.phenotype[i])
            local_2_opt!(route, problem_instance, total_objective) # Improve
            chromosome.phenotype[i] = map(p -> p.id, route)
        end

        # Add chromosome to population
        compute_fitness!(chromosome, problem_instance)
        compute_unfitness!(chromosome, problem_instance)
        chromosome.patient_sets = [Set(route) for route in filter(!isempty, chromosome.phenotype)]
        push!(population, chromosome)
    end

    println("Finished initializing population\n")

    for (i, individual) in enumerate(population)
        individual.id = i
    end

    return population
end


function island_algorithm(n_islands::Int, n_individuals::Int, n_generations::Int, exchange_frequency::Int, n_exchange_individuals::Int, problem_instance::ProblemInstance, mutation_rate::Float64, lns_frequency::Int=10000, verbose::Int=1)

    n_nurses = problem_instance.n_nurses

    # Initialize islands
    islands = Vector{Vector{Chromosome}}()
    @threads for i in 1:n_islands
        push!(islands, initialize_population(n_individuals, n_nurses, problem_instance, 0))
    end

    println("Optimizing ", n_islands, " islands with ", n_individuals, " individuals each on ", nthreads(), " threads\n")

    for i in 1:exchange_frequency:n_generations

        println("----------- Generation ", i - 1, " of ", n_generations, " -----------")
        for j in 1:n_islands
            println("\n--- Island ", j, " ---")
            best_chromosome = sort(islands[j], by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[1]
            println("Top 1 - Fitness: ", @sprintf("%.2f", best_chromosome.fitness), " Time unfitness: ", @sprintf("%.2f", best_chromosome.time_unfitness), " Strain unfitness: ", @sprintf("%.2f", best_chromosome.strain_unfitness))
            println("Number of unique individuals: ", count_unique_individuals(islands[j]), "/", n_individuals)
            println("")
        end
        println("-----------------------------------------\n")


        @threads for j in 1:n_islands
            # Perform genetic algorithm on each island
            islands[j] = genetic_algorithm(islands[j], problem_instance, n_individuals, exchange_frequency, mutation_rate, n_point_crossover, n_nurses, lns_frequency, 0)
        end

        visited = Set{Int}()
        # Exchange individuals between islands
        for j in 1:n_islands
            # Next island
            next = mod1(j + 1, n_islands)
            while next in visited
                next = mod1(next + 1, n_islands)
            end
            # Exchange n_exchange_individuals random individuals
            for k in 1:n_exchange_individuals
                # Remove random individual from current island
                exchange_individual = splice!(islands[j], rand(eachindex(islands[j])))
                # Add to next island
                push!(islands[next], exchange_individual)
            end
            push!(visited, next)
        end

    end

    println("Performing genetic algorithm on merged population")

    # Merge islands
    population = Vector{Chromosome}()
    for island in islands
        sorted_island = sort(island, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))
        population = vcat(population, sorted_island[1:floor(Int, n_individuals / n_islands)+1])
    end
    population = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[1:n_individuals]

    # Perform genetic algorithm on merged population
    population = genetic_algorithm(population, problem_instance, n_individuals, floor(Int, n_generations / 2), mutation_rate, n_point_crossover, n_nurses, verbose)

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

