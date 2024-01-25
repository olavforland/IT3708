

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