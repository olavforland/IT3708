module Objective

export time_unfitness_objective, strain_unfitness_objective, time_fitness_objective, total_objective

using ..DataParser: ProblemInstance, Patient

function total_objective(route::Vector{Patient}, instance::ProblemInstance)
    return (time_unfitness_objective(route, instance), strain_unfitness_objective(route, instance), time_fitness_objective(route, instance))
end

function time_unfitness_objective(route::Vector{Patient}, instance::ProblemInstance)

    elapsed_time = 0.0
    time_violation = 0.0

    prev_patient = 1
    for patient in route
        # Add travel time
        elapsed_time += instance.travel_times[prev_patient][patient.id+1]
        # If arrive early, wait
        elapsed_time += max(patient.start_time - elapsed_time, 0)
        # If arrive late, add to time violation
        time_violation += max(elapsed_time - (patient.end_time - patient.care_time), 0)
        # Add care time
        elapsed_time += patient.care_time
        # Increment by one due to 1-indexing
        prev_patient = patient.id + 1
    end

    return time_violation
end

function strain_unfitness_objective(route::Vector{Patient}, instance::ProblemInstance)

    nurse_strain = 0
    for patient in route
        nurse_strain += patient.demand
    end

    return max(nurse_strain - instance.nurse_capacity, 0)

end

function time_fitness_objective(route::Vector{Patient}, instance::ProblemInstance)

    travel_times = instance.travel_times
    travel_time = 0.0
    # Start at depot
    prev_patient = 1
    for patient in route
        # Increment by one due to 1-indexing
        travel_time += travel_times[prev_patient][patient.id+1]
        # Increment by one due to 1-indexing
        prev_patient = patient.id + 1
    end
    # End at depot
    travel_time += travel_times[prev_patient][1]

    return travel_time

end


end # module
