import numpy as np

import seaborn as sns

from abc import ABC, abstractmethod

class SGA(ABC):

    def __init__(self, size, bits) -> None:
        super().__init__()

        self.size = size
        self.bits = bits
        self.population = np.random.choice([0, 1], size=(size, bits))


    def decode(self, population):

        bit_position_values = 2**np.arange(0, self.bits)[::-1]
        decoded = population @ bit_position_values
        return decoded / 2**(self.bits-7)


    def _roulette_wheel_probas(self, fitness):

        probas = (fitness - fitness.min() + 1e-6) / (fitness - fitness.min() + 1e-6).sum()
        return probas

    def select_parents(self, population, fitness, maximize=True):


        probas = self._roulette_wheel_probas(fitness)

        if not maximize:
            probas = (1 - probas + 1e-6) / (1 - probas + 1e-6).sum()

        parent_ind = np.random.choice(self.size, self.size, replace=True, p=probas)
        return population[parent_ind]


    def crossover(self, parents, p_crossover=0.6):
        size, bits = parents.shape

        crossover_points = np.random.randint(1, bits, size=size)    
        
        new_parents = parents.copy()
        old_parents = parents.copy()

        for i in range(0, size-1, 2):

            if np.random.rand() > p_crossover:
                continue
            
            crossover_point = crossover_points[i]

            new_parents[i, :crossover_point] = old_parents[i+1, :crossover_point]
            new_parents[i+1, :crossover_point] = old_parents[i, :crossover_point]

        return new_parents


    def mutate(self, offspring, p_mut=0.05):

        mutation = np.random.choice([0, 1], size=offspring.shape, p=[1-p_mut, p_mut])
        mutated_offspring = np.logical_xor(offspring, mutation)
        return mutated_offspring
        
        
    def select_survivors(self, parents, offspring, fitness, age_based=True, maximize=True):
        if age_based:
            return offspring
        
        population = np.concatenate([parents, offspring], axis=0)
        if maximize:
            indices = fitness.argsort()[::-1][:self.size]
        else:
            indices = fitness.argsort()[:self.size]

        return population[indices]
    

    @abstractmethod
    def get_fitness(self, population, **kwargs):
        pass
        
    
    def run(self, generations=10, maximize=True, p_crossover=0.6, p_mut=0.05, **kwargs):
        
        mean_fitness = []
        best_fitness = []
        populations = []

        g = 0

        population = self.population.copy()

        
        while g < generations + 1:
            
            fitness = self.get_fitness(population, **kwargs)
            mean_fitness.append(fitness.mean())
            if maximize:
                best_fitness.append(fitness.max())
            else:
                best_fitness.append(fitness.min())

            populations.append(population)
            
            parents = self.select_parents(population, fitness, maximize=maximize)

            np.random.shuffle(parents)

            offspring = self.crossover(parents, p_crossover=p_crossover)
            
            offspring = self.mutate(offspring, p_mut=p_mut)

            fitness = self.get_fitness(np.concatenate([parents, offspring], axis=0), **kwargs)

            population = self.select_survivors(parents, offspring, fitness=fitness, age_based=False, maximize=maximize)

            g += 1

        

        return populations, mean_fitness, best_fitness

class SGA_Sine(SGA):
    
    def __init__(self, size, bits) -> None:
        super().__init__(size, bits)


    def get_fitness(self, population, **kwargs):
        population_values = self.decode(population)

        lb, ub = kwargs.get('feasible_region', (0, 128))
        penalty = 0
        lower_penalty = np.where(
            population_values < lb, 
            lb - population_values, 
            0
        )
        upper_penalty = np.where(
            population_values > ub, 
            population_values - ub,
            0
        )
        penalty = (lower_penalty + upper_penalty)

        return np.sin(population_values) - penalty
    
class SGA_LinReg(SGA):

    def __init__(self, size, bits, linreg) -> None:
        super().__init__(size, bits)

        self.linreg = linreg

    def get_fitness(self, population, **kwargs):
        fitness = np.zeros(population.shape[0])

        for i, feature_mask in enumerate(population):
            fitness[i] = self.linreg.get_fitness(feature_mask, rng=42)

        return fitness