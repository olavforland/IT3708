using Images
using Colors

function read_image(file_path::String)
    # Read image from file path and return a Vector of Vector of RGB values
    # Images are in folder Project3/images
    # path must be relative to the Project3 folder
    # Example: read_image("images/XXXXX/Test image.jpg")


    img = load(file_path)
    rgb_array = Float64.(channelview(img))
    _, height, width = size(rgb_array)
    println("Loaded image of size \n height:", height, "\n width: ", width, "\n total pixels: ", height * width, "\n")

    #convert to vector of vectors of tuples of RGB values
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


#img = read_image("images/86016/Test image.jpg")
#print(img[1:3])


