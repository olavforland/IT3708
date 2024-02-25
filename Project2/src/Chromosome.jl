module Chromosome

struct Chromosome
    genes::Vector{Int}
    fitness::Float64
    unfitness::Float64
end

end # module