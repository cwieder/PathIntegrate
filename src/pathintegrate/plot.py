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

    block_importance = pd.DataFrame(pi_model.A_corrected_, index=pi_model.omics_names, columns=range(0, len(pi_model.A_corrected_)))
    sns.heatmap(
        data=block_importance,
        cmap='Blues',
        annot=True,
        fmt='.2f',
        linewidths=.5,
        square=True)
    plt.ylabel('Omics view')
    plt.xlabel('Percentage explained in Y per LV')
    plt.title('Omics view importance')

    plt.tight_layout()
    if outfile:
        plt.savefig(outfile)
    plt.show()

