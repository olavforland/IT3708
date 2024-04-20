module EvalInterface


using ..Genetics: Chromosome, get_segment_mask
using ..Problem: ProblemInstance


function draw_segments(chromosome::Chromosome, instance::ProblemInstance, path::String="output.txt")
    mask = get_segment_mask(chromosome, instance)
    h = instance.height
    w = instance.width
    img = [[255 for _ in 1:w] for _ in 1:h]

    for i in 1:h-1
        for j in 1:w-1
            if (mask[i][j] != mask[i+1][j]) || (mask[i][j] != mask[i][j+1])
                img[i][j] = 0
            end #if
        end #for
    end #for

    #Draw edges of image
    for i in 1:h
        img[i][1] = 0
        img[i][w] = 0
    end #for
    for j in 1:w
        img[1][j] = 0
        img[h][j] = 0
    end #for

    pre_path = "evaluator/student_segments/"
    #write vector to txt file
    io = open(pre_path * path, "w")
    for line in img
        println(io, join(line, ","))
    end #for
    close(io)
end #function


end #module