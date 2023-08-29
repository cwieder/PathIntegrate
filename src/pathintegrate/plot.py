# plotting functions

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def omics_view_importance(pi_model, outfile=None):
    """Plot the importance of each omics view in the MultiView model
    
    Args:
        pi_model (object): Fitted PathIntegrate multi-view model.
    """
    plt.figure(figsize=(7, 5))
    block_importance = pd.DataFrame(pi_model.A_corrected_*100, index=pi_model.omics_names, columns=range(0, pi_model.A_corrected_.shape[1]))
    pev_per_lv = np.round(pi_model.explained_var_y_, 2)*100
    sns.heatmap(
        data=block_importance,
        cmap='Blues',
        annot=True,
        fmt='.2f',
        linewidths=.5,
        square=True,
        vmin=0,
        vmax=100)
    
    plt.ylabel('Omics view')
    plt.xticks([i+0.5 for i in range(0, pi_model.A_corrected_.shape[1])], [str(n+1) + ': ' + str(i)+'%' for n, i in enumerate(pev_per_lv)])
    plt.xlabel('Percentage explained in Y per LV')
    plt.title('Omics view importance')

    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
    plt.show()


# def top_n_pathways(pi_model, n_top_paths=10, outfile=None):

#     if pi_model.mv:
#         pass
        
#     elif pi_model.sv:
#         paths_df = pi_model.feature_importances_.sort_values(ascending=False).head(n_top_paths)
#         sns.barplot()

#     plt.tight_layout()
#     if outfile:
#         plt.savefig(outfile)    
#     plt.show()
