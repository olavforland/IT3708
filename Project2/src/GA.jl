
module GA

export initialize_population, genetic_algorithm

using Random

using ..DataParser
using ..Genetics
using ..Crossover
using ..Mutation
using ..Selection
using ..TSPHeuristic
using ..VNSHeuristic: construct_solution!, improve_solution!
using ..Utils: count_unique_individuals


function genetic_algorithm(problem_instance::ProblemInstance, n_individuals::Int, n_generations::Int, mutation_rate::Float64, n_nurses::Int)
    population = initialize_population(n_individuals, n_nurses, problem_instance)

    for generation in 1:n_generations
        p1, p2 = tournament_selection(population, 2)
        c1, c2 = two_point_crossover(p1, p2, n_nurses)
        swap_mutation!(c1, mutation_rate), swap_mutation!(c2, mutation_rate)

        # Perform routing step
        construct_solution!(problem_instance, c1)
        construct_solution!(problem_instance, c2)
        
        # improve_solution!(problem_instance, c1)
        # improve_solution!(problem_instance, c2)
        
        # Evaluate fitness and unfitness
        compute_fitness!(c1, problem_instance), compute_fitness!(c2, problem_instance)
        compute_unfitness!(c1, problem_instance), compute_unfitness!(c2, problem_instance)


        # Get favoured offspring  based on fitness
        child = c1.fitness < c2.fitness ? c1 : c2

        # Check if child in population
        
        if generation % 1000 == 0
            best_chromosome = sort(population, by=p -> (p.time_unfitness, p.strain_unfitness, p.fitness))[1]
            println("Generation: ", generation, " Fitness: ", best_chromosome.fitness, " Time unfitness: ", best_chromosome.time_unfitness, " Strain unfitness: ", best_chromosome.strain_unfitness)
            println("Number of unique individuals: ", count_unique_individuals(population), "/", n_individuals)
            println("Number of nurses used: ", length(unique(best_chromosome.genotype)))
            println("")
            
            #println("Child 1: ", c1.fitness, " ", c1.time_unfitness, " ", c1.strain_unfitness)
            # println("Child 2: ", c2.fitness, " ", c2.time_unfitness, " ", c2.strain_unfitness)
        end
        if in(join(child.genotype, ","), map(c -> join(c.genotype, ","), population))
            continue
        end

        survivor_selection!(population, child)

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

    shuffled_patients = shuffle(1:length(problem_instance.patients))
    for i in 1:n_individuals
        delta = 31#shuffled_patients[i]
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

            construct_solution!(problem_instance, chromosome)
            # improve_solution!(problem_instance, chromosome)

            compute_unfitness!(chromosome, problem_instance)
            
            if (chromosome.time_unfitness > prev_time_unfitness || chromosome.strain_unfitness > prev_strain_unfitness)
                current_nurse += 1
                prev_time_unfitness = chromosome.time_unfitness
                prev_strain_unfitness = chromosome.strain_unfitness
            
            end
        end


        # Find most violated routes and distribute half the patients on a new nurse
        routes_most_violated = sort(1:length(chromosome.phenotype), by=p -> (chromosome.route_strain_unfitness[p], chromosome.route_time_unfitness[p]), rev=true)
        for nurse in routes_most_violated
            if chromosome.route_strain_unfitness[nurse] == 0 && chromosome.route_time_unfitness[nurse] == 0
                continue
            end
            current_nurse += 1
            patients_sorted = sort(chromosome.phenotype[nurse], by = p -> ranked_patients[p].rank)
            
            # Assign the first half of the patients to the current nurse
            for patient in patients_sorted[1:div(length(patients_sorted), 2)]
                chromosome.genotype[patient] = current_nurse
            end

            if current_nurse + 1 > n_nurses
                break
            end
            
        end

        # end
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