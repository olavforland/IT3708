module Problem

using ..Utils: knn_dict, read_image

export ProblemInstance

struct ProblemInstance
    pixels::Vector{Vector{Tuple{Float64,Float64,Float64}}}
    knn::Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}
    height::Int
    width::Int

    function ProblemInstance(path::String, n_neighbors::Int=4)
        pixels = read_image(path)
        knn = knn_dict(pixels, n_neighbors)
        height = length(pixels)
        width = length(pixels[1])
        println("Problem instance created")
        new(pixels, knn, height, width)
    end
end

end #module