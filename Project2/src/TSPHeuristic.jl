module TSPHeuristic

export nearest_neighbor_heuristic

using ..DataParser: ProblemInstance, Patient
using ..Genetics: Chromosome




# function to generate a solution based on the nearest neighbor heuristic
function nearest_neighbor_heuristic(instance::ProblemInstance, chromosome::Chromosome, nurse::Int)
    genotype = chromosome.genotype
    nurse_patient_idx = findall(x -> x == nurse, genotype)
    # if nurse_patient_idx == []
    #     return []
    # end
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
            if travel_times[current+1][patient.id+1] < min_time && ((current, patient.id) ∉ instance.inadmissable_presedence)
                min_time = travel_times[current+1][patient.id+1]
                nearest_patient = patient
            end
        end
        if isnothing(nearest_patient)
            # Add the remaining nurse patients to the greedy order
            for patient in nurse_patients
                if travel_times[current+1][patient.id+1] < min_time
                    min_time = travel_times[current+1][patient.id+1]
                    nearest_patient = patient
                end
            end
        end

        # add the nearest patient to the greedy order
        push!(greedy_order, nearest_patient.id)

        # remove the nearest patient from the list of patients left to visit
        deleteat!(nurse_patients, findall(x -> x == nearest_patient, nurse_patients))

        # set the current patient to the nearest patient
        current = nearest_patient.id
    end # end while

    return greedy_order
end # end function nearest_neighbor_heuristic



function savelsbergh_heuristic(instance::ProblemInstance, chromosome::Chromosome, nurse::Int)
    genotype = chromosome.genotype
    nurse_patient_idx = findall(x -> x == nurse, genotype)
    nurse_patients = instance.patients[nurse_patient_idx]
    travel_times = instance.travel_times
    depot_coords = instance.depot_coords
    nurse_capacity = instance.nurse_capacity
    depot_return_time = instance.depot_return_time

    #Making two depot patients 
    depot_start = Patient(0, depot_coords[1], depot_coords[2], 0, 0, 0, 0)
    depot_end = Patient(0, depot_coords[1], depot_coords[2], 0, depot_return_time, depot_return_time, 0)

    route = [depot_start, depot_end]

    while !isempty(nurse_patients)

        optimal_insertion_places = criterion_one_first_phase(nurse_patients, travel_times, depot_coords, nurse_capacity, depot_return_time, route)

        if isnothing(optimal_insertion_places)
            return [p.id for p in route]
            # TODO: TOCONSIDER
            # Not sure if it should break or return current path at this point.
        end

        best_patient, insertion_spots = criterion_two_first_phase(route, travel_times, optimal_insertion_places)
        insert!(route, insertion_spots[2], best_patient)
        deleteat!(nurse_patients, findall(x -> x == best_patient, nurse_patients))

    end #end while

    #remove depot patients
    deleteat!(route, 1)
    deleteat!(route, length(route))
    return [p.id for p in route]
end # end function savelsbergh_heuristic


#----------------- Helper functions -----------------#

function criterion_one_first_phase(unrouted_patients::Vector{Patient}, travel_times::Vector{Vector{Float64}}, depot_coords::Tuple{Int,Int}, nurse_capacity::Int, depot_return_time::Int, current_route::Vector{Patient})
    #Find the best insertion spot for each patient
    # i : previous patient
    # j : next patient
    # u : new patient


    insertion_spots = Dict{Patient,Tuple{Int,Int}}()
    for patient in unrouted_patients
        best_spot = nothing
        h = -Inf
        for i in 1:length(current_route)-1
            if !cfi(current_route, i, i + 1, patient, nurse_capacity, depot_return_time, travel_times)
                continue
            end # end if
            # For all feasible insertion spots
            l_u = patient.end_time - patient.care_time #latest start time for patient
            D_i = departure_time(travel_times, current_route, i) #departure time from previous patient
            t_i_u = travel_times[current_route[i].id+1][patient.id+1] #travel time from previous patient to unrouted patient
            pfs_i = pfs(current_route, i, patient, travel_times) #possible forward shift in time of the departure time at i causing no violation of the time-window constraints along the path (i, ... , i + l);
            emc_ = emc(i, patient, i + 1, travel_times, current_route) #extra milage cost
            e_u = patient.start_time #Earliest start time for patient
            h_patient = min(l_u - max(D_i + t_i_u, e_u), pfs_i - emc_)
            if h_patient > h
                best_spot = i
                h = h_patient
            end # end if
            insertion_spots[patient] = (best_spot, best_spot + 1)
        end # end for i
    end # end for patient

    return insertion_spots
end # end function criterion_one_first_phase

function criterion_two_first_phase(current_route::Vector{Patient}, travel_times::Vector{Vector{Float64}}, insertion_spots::Dict{Patient,Tuple{Int,Int}})
    #Find the best patient to insert in the route
    best_patient = nothing
    h = Inf
    for (p, value) in pairs(insertion_spots)
        h_insertion = emc(value[1], p, value[2], travel_times, current_route)
        if h_insertion < h
            best_patient = p
            h = h_insertion
        end # end if
    end # end for
    return best_patient, insertion_spots[best_patient]
end # end function criterion_two_first_phase


#possible forward shift
#possible forward shift in time of the departure time at j causing no violation of the time-window constraints along the path (j, ... , i + l);
function pfs(current_route::Vector{Patient}, i::Int, patient::Patient, travel_times::Vector{Vector{Float64}})
    possible_new_route = copy(current_route)
    insert!(possible_new_route, i + 1, patient)
    k = i
    D_k = departure_time(travel_times, possible_new_route, k) #departure time from patient k, need to check which patient k is
    l_N = last(current_route).end_time
    if length(current_route) == 2
        SIGMA = 0
    else
        SIGMA = sum(p.care_time for (index, p) in enumerate(current_route) if index >= k && index < length(current_route)) + sum(travel_times[p.id+1][current_route[index+1].id+1] for (index, p) in enumerate(current_route) if index >= k && index < length(current_route))
    end

    return l_N - D_k - SIGMA
end # end function pfs

# extra milage cost
function emc(i::Int, u::Patient, j::Int, travel_times::Vector{Vector{Float64}}, current_route::Vector{Patient})
    # i : previous patient
    # j : next patient
    # u : new patient
    new = u
    prev = current_route[i]
    next = current_route[j]
    possible_new_route = copy(current_route)
    insert!(possible_new_route, i + 1, u)
    D_i = departure_time(travel_times, current_route, i) #departure time from previous patient
    t_i_u = travel_times[prev.id+1][new.id+1] #travel time from previous patient to unrouted patient
    t_u_j = travel_times[new.id+1][next.id+1] #travel time from unrouted patient to next patient
    A_j = departure_time(travel_times, possible_new_route, i + 1) + travel_times[new.id+1][next.id+1] #arrival time at next patient
    return max(D_i + t_i_u, new.start_time) + t_u_j - A_j
end # end function emc

#check feasibility of insertion
function cfi(current_route::Vector{Patient}, i::Int, j::Int, u::Patient, nurse_capacity::Int, depot_return_time::Int, travel_times::Vector{Vector{Float64}})
    # i : previous patient
    # j : next patient
    # u : new patient
    if sum([i.demand for i in current_route]) + u.demand > nurse_capacity
        return false
    end # end if

    possible_new_route = copy(current_route)
    insert!(possible_new_route, i + 1, u)


    elapsed_time = 0.0
    time_violation = 0.0
    prev_patient = 1
    for patient in possible_new_route
        # Add travel time
        elapsed_time += travel_times[prev_patient][patient.id+1]
        # If arrive early, wait
        elapsed_time += max(patient.start_time - elapsed_time, 0)
        # If arrive late, add to time violation
        time_violation += max(elapsed_time - (patient.end_time - patient.care_time), 0)
        # Add care time
        elapsed_time += patient.care_time
        # Increment by one due to 1-indexing
        prev_patient = patient.id + 1
    end
    if time_violation > 0
        return false
    end # end if

    return true
end # end function cfi

function departure_time(travel_times::Vector{Vector{Float64}}, current_route::Vector{Patient}, i::Int)
    # departure time from patient i
    time = 0
    for k in 1:i-1
        time += travel_times[current_route[k].id+1][current_route[k+1].id+1]
        if current_route[k].start_time > time
            time = current_route[k].start_time
        end # end if
        time += current_route[k].care_time
    end # end for
    return time
end # end function departure_time


end # module TSPHeuristic


