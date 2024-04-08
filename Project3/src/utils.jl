using Images

function read_image(file_path)::Vector{Vector{Tuple{Int,Int,Int}}}
    # Read image from file path and return a 2D array of RGB values
    # Images are in folder Project3/images
    # path must be relative to the Project3 folder
    # Example: read_image("images/XXXXX/Test image.jpg")


    img = load(file_path)
    rgb_vals = [[(r, g, b) for (r, g, b) in eachrow(row)] for row in eachslice(img)]
    return rgb_vals
end


