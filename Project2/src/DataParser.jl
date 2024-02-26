module DataParser

export parse_data, Patient, ProblemInstance

using JSON

mutable struct Patient
    id::Int
    x_coord::Int
    y_coord::Int    
    demand::Int
    start_time::Int
    end_time::Int
    care_time::Int
    angle_from_depot::Union{Float64, Nothing}
    rank::Union{Int, Nothing}

    # Constructor that sets angle_from_depot and rank to nothing by default
    Patient(id::Int, x_coord::Int, y_coord::Int, demand::Int, start_time::Int, end_time::Int, care_time::Int) = new(id, x_coord, y_coord, demand, start_time, end_time, care_time, nothing, nothing)
end

struct ProblemInstance
    depot_return_time::Int
    depot_coords::Tuple{Int, Int}
    n_nurses::Int
    nurse_capacity::Int
    patients::Vector{Patient}
    travel_times::Vector{Vector{Float64}}
end

function parse_data(filepath::String)
    # Read JSON file
    raw_data = JSON.parsefile(filepath)

    # Extract data for depot
    depot_data = raw_data["depot"]
    depot_return_time = depot_data["return_time"]
    depot_coords = (depot_data["x_coord"], depot_data["y_coord"])

    # Extract data for patients
    patients_data = raw_data["patients"]
    patients = [Patient(parse(Int, id), pat["x_coord"], pat["y_coord"], pat["demand"],
                        pat["start_time"], pat["end_time"], pat["care_time"]) for (id, pat) in patients_data]

    rank_patients!(patients, depot_coords) 

    travel_times_data = raw_data["travel_times"]
    travel_times = Vector{Vector{Float64}}()
    for row in travel_times_data
        push!(travel_times, [Float64(elem) for elem in row])
    end

    n_nurses = raw_data["nbr_nurses"]
    nurse_capacity = raw_data["capacity_nurse"]

    # Return all parsed data
    return ProblemInstance(depot_return_time, depot_coords, n_nurses, nurse_capacity, patients, travel_times)
end


# ------------- Helpers ------------- #


function calculate_angle(patient::Patient, depot::Tuple{Int, Int})
    delta_x = patient.x_coord - depot[1]
    delta_y = patient.y_coord - depot[2]
    raw_angle = atan(delta_x, delta_y)

    # First quadrant
    if delta_x >= 0 && delta_y >= 0
        return pi / 2 - raw_angle
    
    # Second or third quadrant
    elseif delta_x < 0
        return pi / 2 - raw_angle
    end

    # Fourth quadrant
    return 5 * pi / 2 - raw_angle
end

function rank_patients!(patients::Vector{Patient}, depot::Tuple{Int, Int})
    for patient in patients
        patient.angle_from_depot = calculate_angle(patient, depot)
    end

    sort!(patients, by = patient -> patient.angle_from_depot)

    for (i, patient) in enumerate(patients)
        patient.rank = i
    end
end

end # module