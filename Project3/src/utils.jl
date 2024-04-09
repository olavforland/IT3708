module Utils

using Images
using Colors
using DataStructures

export read_image, euclidean_distance, min_spanning_tree, mst_to_genotype


function read_image(file_path::String)::Vector{Vector{Tuple{Float64,Float64,Float64}}}

    """
    Read image from file path and return a Vector of Vector of RGB values
    Images are in folder Project3/images
    path must be relative to the Project3 folder
    Example: read_image("images/XXXXX/Test image.jpg")
    """

    img = load(file_path)
    rgb_array = Float64.(channelview(img))
    _, height, width = size(rgb_array)
    num_pixels = height * width
    println("Loaded image of size \n height:", height, "\n width: ", width, "\n total pixels: ", num_pixels, "\n")


    pixel_vectors = Vector{Vector{Tuple{Float64,Float64,Float64}}}()
    for i in 1:height
        row = Vector{Tuple{Float64,Float64,Float64}}()
        for j in 1:width
            push!(row, (rgb_array[1, i, j], rgb_array[2, i, j], rgb_array[3, i, j]))
        end #for
        push!(pixel_vectors, row)
    end #for



    return pixel_vectors
end #read_image


function euclidean_distance(p1::Tuple{Float64,Float64,Float64}, p2::Tuple{Float64,Float64,Float64})::Float64
    """
    function that calculates the euclidean distance between two pixels in rgb-space

    args: two pixels in rgb-space
    returns: float, euclidean distance between the two pixels
    """
    eucl = sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2 + (p1[3] - p2[3])^2)
    return eucl
end #euclidean_distance



function min_spanning_tree(rgbs::Vector{Vector{Tuple{Float64,Float64,Float64}}})::Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}
    """
    function that creates a minimum spanning tree over the image 
    where the nodes are pixels and the edges are the difference 
    in RGB values. The weight of an edge is the euclidean distance.
    The root of the tree is selected randomly to generate a different 
    tree each time.

    args: Vector of Vector of RGB values
    returns: Dict of Tuple of Int, Int to Set of Tuple of Int, Int
    
    """

    h = length(rgbs)
    w = length(rgbs[1])

    num_nodes = h * w

    connected = Set{Tuple{Int,Int}}()

    mst = Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}()
    pq = PriorityQueue{Tuple{Tuple{Int,Int},Tuple{Int,Int}},Float64}()

    is_valid_index = (i, j) -> 1 <= i <= h && 1 <= j <= w

    #randomly select root
    root = (rand(1:h), rand(1:w))
    print("root: ", root, "\n")
    enqueue!(pq, (root, root) => 0.0)

    while !isempty(pq)
        (u, v) = dequeue!(pq)
        if !(v in connected)
            push!(connected, v)
            if !(u in keys(mst))
                mst[u] = Set{Tuple{Int,Int}}()
            end #if
            push!(mst[u], v)
            for (di, dj) in [(0, 1), (1, 0), (0, -1), (-1, 0)]
                if is_valid_index(v[1] + di, v[2] + dj)
                    neighbor = (v[1] + di, v[2] + dj)
                    if !(neighbor in connected)
                        weight = euclidean_distance(rgbs[v[1]][v[2]], rgbs[neighbor[1]][neighbor[2]])
                        enqueue!(pq, (v, neighbor) => weight)
                    end #if
                end #if
            end #for
        end #if
    end #while


    println("Number of nodes in the minimum spanning tree: ", length(connected), "\n")



    return mst


end #min_spanning_tree




function mst_to_genotype(mst::Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}, dims::Tuple{Int,Int})::Vector{Char}
    """
    function that converts a minimum spanning tree to a genotype
    where the edges are represented by characters. 
    The characters are u, d, l, r, n for up, down, left, right, none.

    In for this to work as is, mst has to represent a directed graph
    where the keys are the parent nodes and the values are the children nodes.

    args: Dict of Tuple of Int, Int to Set of Tuple of Int, Int
    returns: Vector of Char
    """

    #initialize genotype as Vector of Vector of Char with dimensions dims
    h, w = dims
    genotype = Vector{Vector{Char}}([fill('n', w) for _ in 1:h])

    for (r, neighbors) in mst
        for v in neighbors
            if r[1] == v[1] + 1
                genotype[v[1]][v[2]] = 'u'
            elseif r[1] == v[1] - 1
                genotype[v[1]][v[2]] = 'd'
            elseif r[2] == v[2] + 1
                genotype[v[1]][v[2]] = 'l'
            elseif r[2] == v[2] - 1
                genotype[v[1]][v[2]] = 'r'
            end #if
        end #for
    end #for

    #flatten genotype
    return return reduce(vcat, genotype)

end #mst_to_genotype


end #module



