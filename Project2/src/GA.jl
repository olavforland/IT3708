
module GA

export initialize_population, genetic_algorithm

using ..DataParser
using ..Genetics
using ..Crossover
using ..Mutation
using ..Selection


function genetic_algorithm(problem_instance::ProblemInstance, n_individuals::Int, n_generations::Int, mutation_rate::Float64, n_nurses::Int)
    population = initialize_population(n_individuals, n_nurses, problem_instance)

    for generation in 1:n_generations
        p1, p2 = tournament_selection(population, 2)
        c1, c2 = two_point_crossover(p1, p2)
        swap_mutation!(c1, mutation_rate), swap_mutation!(c2, mutation_rate)

        # Perform routing step
        c1.phenotype = simple_routing_sort_start_times(c1, problem_instance.patients, n_nurses)
        c2.phenotype = simple_routing_sort_start_times(c2, problem_instance.patients, n_nurses)
        
        # Evaluate fitness and unfitness
        compute_fitness!(c1, problem_instance), compute_fitness!(c2, problem_instance)
        compute_unfitness!(c1, problem_instance), compute_unfitness!(c2, problem_instance)

        # Get favoured offspring  based on fitness
        child = c1.fitness > c2.fitness ? c1 : c2

        # Check if child in population
        if in(child, population)
            continue
        end

        survivor_selection!(population, child) 
        
        if generation % 50000 == 0
            best_chromosome = population[findmin(getfield.(population, :fitness))[2]]
            println("Generation: ", generation, " Fitness: ", best_chromosome.fitness, " Time unfitness: ", best_chromosome.time_unfitness, " Strain unfitness: ", best_chromosome.strain_unfitness)
        end

    end

    return population
end


function initialize_population(n_individuals::Int, n_nurses::Int, problem_instance::ProblemInstance)
    population = Vector{Chromosome}()

    n_patients = length(problem_instance.patients)

    patients = sort(problem_instance.patients, by = p -> p.id)
    ranked_patients = sort(problem_instance.patients, by = p -> p.rank)

    # Initial nurse areas
    nurse_areas = create_initial_nurse_areas(n_nurses, n_patients)
    for i in 1:n_individuals
        delta = rand(1:length(problem_instance.patients))
        chromosome = Chromosome([0 for _ in 1:n_patients])
        chromosome.phenotype = [Vector{Int}() for _ in 1:n_nurses]
        current_nurse = 1

        prev_time_unfitness = 0.0
        prev_strain_unfitness = 0.0

        for j in 1:n_patients
            
            # Current patient
            p = mod1(j + delta, n_patients)
            patient = ranked_patients[p]

            # Add patient to the genotype
            chromosome.genotype[patient.id] = current_nurse

            # Sort the patients by start of time window
            chromosome.phenotype[current_nurse] = push!(chromosome.phenotype[current_nurse], patient.id)
            chromosome.phenotype[current_nurse] = sort(chromosome.phenotype[current_nurse], by = p -> patients[p].start_time)

            compute_unfitness!(chromosome, problem_instance)

            if chromosome.time_unfitness > prev_time_unfitness || chromosome.strain_unfitness > prev_strain_unfitness
                current_nurse += 1
                prev_time_unfitness = chromosome.time_unfitness
                prev_strain_unfitness = chromosome.strain_unfitness


            end

            if current_nurse > n_nurses or 
                break
            end



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