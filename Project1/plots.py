
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def plot_fitness(mean_fitness, best_fitness):

    # Initialize the matplotlib figure
    fig = plt.figure(figsize=(20, 8))

    generations = list(range(len(mean_fitness)))

    sns.set_palette('deep')

    sns.lineplot(x=generations, y=mean_fitness, linewidth=4, alpha=0.8, label='Mean Fitness')
    sns.lineplot(x=generations, y=best_fitness, linewidth=4, alpha=0.8, label='Best Fitness')


    plt.ylabel('Fitness', fontsize=22, fontfamily='serif', labelpad=15)
    plt.xlabel('Generations', fontsize=22, fontfamily='serif', labelpad=15)
    plt.xticks(fontsize=15, fontfamily='serif')
    plt.yticks(fontsize=13, fontfamily='serif')

    # plt.gca().set_xscale('log')

    # Increase the font size of the legend
    legend = plt.legend(
        fontsize=14, title='', 
        frameon=False, prop={'family':'serif', 'size':20}, 
        ncol=1, title_fontproperties=dict(family='serif', size=22, weight='bold')
    )
    legend.set_bbox_to_anchor((0.8, 0.7))


    plt.show()


def plot_entropy(entropy_dict):
    fig = plt.figure(figsize=(20, 8))
    colors = sns.color_palette('deep', n_colors=len(entropy_dict))

    for i, (key, value) in enumerate(entropy_dict.items()):
        sns.lineplot(x=list(range(len(value))), y=value, linewidth=4, alpha=0.8, label=key, color=colors[i])
    
    plt.ylabel('Entropy', fontsize=22, fontfamily='serif', labelpad=15)
    plt.xlabel('Generations', fontsize=22, fontfamily='serif', labelpad=15)
    plt.xticks(fontsize=15, fontfamily='serif')
    plt.yticks(fontsize=13, fontfamily='serif')

    legend = plt.legend(
        fontsize=14, title='',
        frameon=False, prop={'family':'serif', 'size':20},
        ncol=1, title_fontproperties=dict(family='serif', size=22, weight='bold')
    )
    legend.set_bbox_to_anchor((0.8, 0.7))

    plt.show()

def compare_fitness(mean_fitness_dict, best_fitness_dict):
    fig = plt.figure(figsize=(20, 8))
    colors = sns.color_palette('deep', n_colors=len(mean_fitness_dict))

    for i, (key, value) in enumerate(mean_fitness_dict.items()):
        sns.lineplot(x=list(range(len(value))), y=value, linewidth=4, alpha=0.8, label=key, color=colors[i])
    
    for i, (key, value) in enumerate(best_fitness_dict.items()):
        sns.lineplot(x=list(range(len(value))), y=value, linewidth=4, alpha=0.8, label=key, color=colors[i], linestyle='--')

    plt.ylabel('Fitness', fontsize=22, fontfamily='serif', labelpad=15)
    plt.xlabel('Generations', fontsize=22, fontfamily='serif', labelpad=15)
    plt.xticks(fontsize=15, fontfamily='serif')
    plt.yticks(fontsize=13, fontfamily='serif')

    legend = plt.legend(
        fontsize=14, title='',
        frameon=False, prop={'family':'serif', 'size':20},
        ncol=1, title_fontproperties=dict(family='serif', size=22, weight='bold')
    )
    legend.set_bbox_to_anchor((0.8, 0.7))

    plt.show()

def plot_sine(population, generation=0):

    fig = plt.figure(figsize=(20, 8))

    sns.set_palette('deep')
    x = np.linspace(0, 128, 10000)
    y = np.sin(x)

    sns.scatterplot(x=population[generation], y=np.sin(population[generation]), color='red', label='Individual', s=200, alpha=1)
    sns.lineplot(x=x, y=y, linewidth=4, alpha=0.5, label='Sine function')

    plt.ylabel('Fitness', fontsize=22, fontfamily='serif', labelpad=15)
    plt.xlabel('Individuals', fontsize=22, fontfamily='serif', labelpad=15)
    plt.xticks(fontsize=15, fontfamily='serif')
    plt.yticks(fontsize=13, fontfamily='serif')
    

    legend = plt.legend(
        fontsize=14, title='', 
        frameon=False, prop={'family':'serif', 'size':20}, 
        ncol=1, title_fontproperties=dict(family='serif', size=22, weight='bold')
    )
    legend.set_bbox_to_anchor((1, 1))

    plt.title(f'Sine Function and Population Individuals after {generation} Generations', fontsize=22, fontfamily='serif', pad=15)