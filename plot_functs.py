import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from cmcrameri import cm

sns.set_style('ticks')

def plot_block_importance(model):
    plt.figure(figsize=(7, 5))
    sns.heatmap(model.A_corrected_*100, cmap='cmc.devon_r', linewidth=5, square=True,
    cbar_kws={'label': '% variance explained in Y'},
    annot=True,
    vmin=0,
    vmax=100)
    plt.yticks([i + 0.5 for i in range(0, len(model.A_corrected_))], model.omics_names)
    plt.xticks([i+0.5 for i in range(0, 5)], [str(n+1) + ': ' + str(np.around(i*100, 1))+'%' for n, i in enumerate(model.explained_var_y_)])
    # plt.xticks([i+0.5 for i in range(0, 5)], range(1, 6))
    plt.xlabel('Component % explained variance in Y')
    plt.ylabel('Omics view')
    plt.tight_layout()
    # plt.savefig('../Figures/MBPLS_COPD_Model1_3o_block_importance.png', dpi=300)
    plt.show()