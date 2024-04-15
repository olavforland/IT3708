module Problem

using .Utils: knn_dict, read_image

export ProblemInstance

struct ProblemInstance
    pixels::Vector{Vector{Tuple{Float64,Float64,Float64}}}
    knn::Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}
    height::Int = length(pixels)
    width::Int = length(pixels[1])

    function ProblemInstance(path::String, n_neighbors::Int=4)
        pixels = read_image(path)
        knn = knn_dict(pixels, n_neighbors)
        new(pixels, knn)
    end
end

end #module