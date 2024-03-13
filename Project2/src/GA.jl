
module GA

export initialize_population, genetic_algorithm

using Random
using Printf
using Base.Threads
using Distributed

using ..DataParser
using ..Genetics
using ..Crossover
using ..Mutation
using ..Selection: tournament_selection, survivor_selection!
using ..TSPHeuristic
using ..VNSHeuristic: construct_solution!, improve_solution!, local_2_opt!, improve_single_route
using ..Utils: count_unique_individuals
using ..LambdaInterchange: lambda_shift_operation, lambda_interchange_operation
using ..Objective
using ..LargeNeighborhoodSearch: tsp_all_routes!



function genetic_algorithm(problem_instance::ProblemInstance, n_individuals::Int, n_generations::Int, mutation_rate::Float64, n_nurses::Int)
    start_time = time()  # Capture start time    

    # Direct initialization
    t0 = time()  # Capture start time
    println("\nInitializing population with " * string(nthreads()) * " threads\n")
    population = initialize_population(n_individuals, n_nurses, problem_instance)
    initialization_time = time() - t0  # Calculate elapsed time


    improvement_time = 0.0
    lambda_shift_time = 0.0
    lambda_interchange_time = 0.0

    similarity_matrix = 

    for generation in 1:n_generations
        p1, p2 = tournament_selection(population, 2)
        c1, c2 = two_point_crossover(p1, p2, problem_instance)
        # c1, c2 = visma_crossover(p1, p2, n_nurses, problem_instance)
        swap_mutation!(c1, mutation_rate), swap_mutation!(c2, mutation_rate)
        # insert_mutation!(c1, mutation_rate), insert_mutation!(c2, mutation_rate)

        # Perform routing step
        construct_solution!(problem_instance, c1, n_nurses)
        construct_solution!(problem_instance, c2, n_nurses)


        # tsp_all_routes!(c1, problem_instance)
        # c1 = lambda_shift_operation(c1, problem_instance)
        # c1 = lambda_interchange_operation(c1, problem_instance)
        for j in 1:length(c1.phenotype)
            route = map(p -> problem_instance.patients[p], c1.phenotype[j])
            local_2_opt!(route, problem_instance, total_objective) # Improve
            c1.phenotype[j] = map(p -> p.id, route)
        end

        # tsp_all_routes!(c2, problem_instance)

        # c2 = lambda_shift_operation(c2, problem_instance)
        # c2 = lambda_interchange_operation(c2, problem_instance)

        for j in 1:length(c2.phenotype)
            route = map(p -> problem_instance.patients[p], c2.phenotype[j])
            local_2_opt!(route, problem_instance, total_objective) # Improve
            c2.phenotype[j] = map(p -> p.id, route)
        end


        # Evaluate fitness and unfitness
        compute_fitness!(c1, problem_instance), compute_fitness!(c2, problem_instance)
        compute_unfitness!(c1, problem_instance), compute_unfitness!(c2, problem_instance)

        # Get favoured offspring  based on fitness
        child = c1.fitness < c2.fitness ? c1 : c2
        # improvement_time += @elapsed begin
        #     improve_solution!(problem_instance, child)
        # end

        if generation % 1000 == 0
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



        if generation % 10000 == 0
            println("\nPerforming large neighborhood search\n")

            @threads for i in eachindex(population)
                # Perform large neighborhood TSP
                tsp_all_routes!(population[i], problem_instance)

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

        end
        # Check if in population    
        if join(child.genotype, ",") ∉ map(c -> join(c.genotype, ","), population)
            survivor_selection!(population, child)
        end

    end


    for chromosome in population
        for i in 1:length(chromosome.phenotype)
            route = map(p -> problem_instance.patients[p], chromosome.phenotype[i])
            chromosome.phenotype[i] = map(p -> p.id, route)
        end
    end
    total_time = time() - start_time
    println("\nInitialization time: ", @sprintf("%.2f", initialization_time), " s (", @sprintf("%.2f", initialization_time / total_time * 100), "% of total time)")
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

    shuffled_patients = shuffle(1:length(patients))

    @threads for i in 1:n_individuals
        delta = shuffled_patients[i]
        current_nurse = nurse_areas[delta]
        # println("current_nurse: ", current_nurse)
        chromosome = Chromosome([0 for _ in 1:n_patients], n_nurses)
        chromosome.phenotype = [Vector{Int}() for _ in 1:n_nurses]

        chromosome.id = i

        for j in 1:n_patients

            # Current patient
            p = mod1(j + delta, n_patients)
            patient = ranked_patients[p]

            chromosome.genotype[patient.id] = current_nurse

            construct_solution!(problem_instance, chromosome, n_nurses)

            compute_unfitness!(chromosome, problem_instance)
            compute_fitness!(chromosome, problem_instance)
            # improve_solution!(problem_instance, chromosome)


            if (chromosome.time_unfitness > 0.0 || chromosome.strain_unfitness > 0.0) #&& (length(chromosome.phenotype[current_nurse]) > n_patients / n_nurses)


                current_nurse += 1
                current_nurse = mod1(current_nurse, n_nurses)

                # If violation not too bad, we accept the infeasible solution into the population
                if chromosome.time_unfitness / chromosome.route_fitness[mod1(current_nurse - 1, n_nurses)] < 0.15 && chromosome.strain_unfitness / patient.demand < 0.15
                    continue
                end
                # If not, we undo the change making the individual infeasible
                chromosome.genotype[patient.id] = current_nurse
                construct_solution!(problem_instance, chromosome, n_nurses)

                # Recompute (un)fitness
                compute_unfitness!(chromosome, problem_instance)
                compute_fitness!(chromosome, problem_instance)

                # prev_time_unfitness = chromosome.time_unfitness
                # prev_strain_unfitness = chromosome.strain_unfitness

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

