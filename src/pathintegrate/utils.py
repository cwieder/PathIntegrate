import pandas as pd
import numpy as np
import pkg_resources
import scipy.stats as stats
import statsmodels.api as sm

def load_example_data(omicstype):
    """
    Loads example datasets

    Args:
        omicstype (str): type of omics for example data. 
            Available options are "metabolomics" or "proteomics". 
            Data are from Su et al 2020 https://doi.org/10.1016/j.cell.2020.10.037.

    Returns:
        pre-processed omics data matrix consisting of m samples and n entities (metabolites/genes) in the form of a pandas DataFrame. 
        Contains one of more metadata columns at the end.
    """

    if omicstype == "metabolomics":
        stream = pkg_resources.resource_stream(__name__, 'data/metabolomics_example.csv')
        f = pd.read_csv(stream, index_col=0, encoding='latin-1')
        return f
    if omicstype == "proteomics":
        stream = pkg_resources.resource_stream(__name__, 'data/proteomics_example.csv')
        f = pd.read_csv(stream, index_col=0, encoding='latin-1')
        return f
