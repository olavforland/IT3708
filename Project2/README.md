

### Requirements
- File parser. Should return a structure containing information about problem:
  - Distance matrix
  - Patient coords (+ polar degrees from depot, maintain mapping from original patient number to new priority)
  - Nurse capacity
- Create Chromosome struct.
  - 1D vector of size equal to number of patients containing integers
  - Needs fitness/unfitness attribute
- Create initial population. Use sweeping based on polar degrees.
- Variations operators. Look into how it is most conveniently handled in Julia.
  - Mutation: Random swap of customers between nurses
  - Crossover: 2-point as a start
- Function that solves TSP with heuristics
  - Start with greedy based on either time or euclidean distance
  - Add local / tabu search
  - Also need exact solver for the final population
- Fitness / Unfitness functions that act on Chromosome struct:
  - Fitness is total time spent travelling
  - Unfitness function maps the violated constraints to a number
  - Maintain cache for seen Chromosomes
- Selection
  - Parent selection (2 x binary tournaments based on (un)fitness)
  - Survivor selection according to the 4 sets introduced in article
