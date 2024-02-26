module TSPHeuristic

using ..DataParser: ProblemInstance, Patient
using ..Genetics: Chromosome


# function to generate a solution based on the nearest neighbor heuristic
function nearest_neighbor_heuristic(instance::ProblemInstance, chromosome::Chromosome, nurse::Int)
    genotype = chromosome.genotype
    nurse_patient_idx = findall(x -> x == nurse, genotype)
    nurse_patients = instance.patients[nurse_patient_idx]
    travel_times = instance.travel_times

    greedy_order = []

    current = 0 # start at the depot, current takes the id of the patient we are currently visiting

    # while there are patients left to visit
    while length(nurse_patients) > 0
        # find the nearest patient
        min_time = Inf
        nearest_patient = nothing
        for patient in nurse_patients
            if travel_times[current+1][patient.id+1] < min_time
                min_time = travel_times[current +1][patient.id +1]
                nearest_patient = patient
            end
        end

        # add the nearest patient to the greedy order
        push!(greedy_order, nearest_patient.id)

        # remove the nearest patient from the list of patients left to visit
        deleteat!(nurse_patients, findall(x -> x == nearest_patient, nurse_patients))

        # set the current patient to the nearest patient
        current = nearest_patient.id
    end

    return greedy_order
end

end # module TSPHeuristic