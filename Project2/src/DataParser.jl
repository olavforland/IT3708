module DataParser

export parse_data

using JSON

struct Patient
    id::Int
    x_coord::Int
    y_coord::Int
    demand::Int
    start_time::Int
    end_time::Int
    care_time::Int
end

struct ProblemInstance
    depot_return_time::Int
    depot_coords::Tuple{Int, Int}
    patients::Vector{Patient}
    travel_times::Vector{Vector{Float64}}
end

# Assume the JSON file is in the data directory and named patients_data.json
function parse_data(filepath::String)
    # Read the JSON file
    raw_data = JSON.parsefile(filepath)

    # Extract data for the depot
    depot_data = raw_data["depot"]
    depot_return_time = depot_data["return_time"]
    depot_coords = (depot_data["x_coord"], depot_data["y_coord"])

    # Extract data for the patients
    patients_data = raw_data["patients"]
    patients = [Patient(parse(Int, id), pat["x_coord"], pat["y_coord"], pat["demand"],
                        pat["start_time"], pat["end_time"], pat["care_time"]) for (id, pat) in patients_data]

    travel_times_data = raw_data["travel_times"]
    travel_times = Vector{Vector{Float64}}()
    for row in travel_times_data
        push!(travel_times, [Float64(elem) for elem in row])
    end

    # Return all parsed data
    return ProblemInstance(depot_return_time, depot_coords, patients, travel_times)
end

end # module