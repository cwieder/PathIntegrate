import pandas as pd
import numpy as np
import sspa
import sklearn
from mbpls.mbpls import MBPLS
import plot_functs
from app import launch_network_app
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA


def sspa_svd(mat, pathway_df, min_entity=2, return_molecular_importance=False):

    """
    Tomfohr et al 2004 SVD/PLAGE method for single sample pathway analysis

    Args:
        mat (pd.DataFrame): pandas DataFrame omics data matrix consisting of m rows (samples) and n columns (entities).
        Do not include metadata columns
        pathways (pd.DataFrame): Dictionary of pathway identifiers (keys) and corresponding list of pathway entities (values).
        Entity identifiers must match those in the matrix columns
        min_entity (int): minimum number of metabolites mapping to pathways for ssPA to be performed
        return_molecular_importance (bool): if True, return molecular importances based on loadings from PCA

    Returns:
        pandas DataFrame of pathway scores derived using the PLAGE method. Columns represent pathways and rows represent samples.
    """
    pathways = sspa.utils.pathwaydf_to_dict(pathway_df)

    pathway_activities = []
    molecular_importance = {}

    # Create pathway matrices
    pathway_ids = []
    for pathway, compounds in pathways.items():
        single_pathway_matrix = mat.drop(mat.columns.difference(compounds), axis=1)

        if single_pathway_matrix.shape[1] >= min_entity:
            pathway_ids.append(pathway)
            pathway_mat = single_pathway_matrix.T.values

            pca = PCA(n_components=1, random_state=0)
            
            pathway_activities.append(pca.fit_transform(single_pathway_matrix)[:, 0])
            loadings = pca.components_[0]
            molecular_importance[pathway] = pd.DataFrame(loadings, index=single_pathway_matrix.columns, columns=['loadings'])

    pathway_activities_df = pd.DataFrame(pathway_activities, columns=mat.index, index=pathway_ids).T
    
    if return_molecular_importance:
        return [pathway_activities_df, molecular_importance]
    else:
        return pathway_activities_df


class PathIntegrate:

    def __init__(self, omics_data:dict, metadata, pathway_source, sspa_scoring='svd', min_coverage=5):
        self.omics_data = omics_data
        self.metadata = metadata
        self.pathway_source = pathway_source
        self.pathway_dict = sspa.utils.pathwaydf_to_dict(pathway_source)
        self.sspa_scoring = sspa_scoring
        self.min_coverage = min_coverage
        
            # CHANGE BACK TO NORMAL SSPA 
        sspa_methods = {'svd': sspa_svd, 'ssGSEA': sspa.sspa_ssGSEA, 'kpca': sspa.sspa_kpca, 'ssClustPA': sspa.sspa_ssClustPA, 'zscore': sspa.sspa_zscore}
        self.sspa_method = sspa_methods[self.sspa_scoring]
        self.sspa_scores_mv = None
        self.sspa_scores_sv = None
        self.coverage = self.get_multi_omics_coverage()

        self.mv = None
        self.sv = None

        self.labels = pd.factorize(self.metadata)[0]
    
    def get_multi_omics_coverage(self):
        all_molecules = sum([i.columns.tolist() for i in self.omics_data.values()], [])
        coverage = {k: len(set(all_molecules).intersection(set(v))) for k, v in self.pathway_dict.items()}
        return coverage

    def MultiView(self, ncomp=2):
        sspa_scores = [self.sspa_method(i, self.pathway_source, self.min_coverage, return_molecular_importance=True) for i in self.omics_data.values()]
        self.sspa_scores_mv = dict(zip(self.omics_data.keys(), [i[0] for i in sspa_scores]))
        # self.sspa_scores_mv = {k: self.sspa_method(v, self.pathway_source, self.min_coverage) for k, v in self.omics_data.items()}
        mv = MBPLS(n_components=ncomp)
        mv.fit([i.copy(deep=True) for i in self.sspa_scores_mv.values()], self.labels)

        # compute VIP and scale VIP across omics
        vip_scores = VIP_multiBlock(mv.W_, mv.Ts_, mv.P_, mv.V_)
        vip_df = pd.DataFrame(vip_scores, index=sum([i.columns.tolist() for i in self.sspa_scores_mv.values()], []))
        vip_df['Name'] = vip_df.index.map(dict(zip(self.pathway_source.index, self.pathway_source['Pathway_name'])))
        vip_df['Source'] = sum([[k] * v.shape[1] for k, v in self.sspa_scores_mv.items()], [])
        vip_df['VIP_scaled'] = vip_df.groupby('Source')[0].transform(lambda x: StandardScaler().fit_transform(x.values[:,np.newaxis]).ravel())

        mv.name = 'MultiView'
        mv.molecular_importances = dict(zip(self.omics_data.keys(), [i[1] for i in sspa_scores]))
        mv.beta = mv.beta_.flatten()
        mv.vip = vip_df
        mv.omics_names = list(self.omics_data.keys())
        mv.sspa_scores = self.sspa_scores_mv
        mv.coverage = self.coverage
        self.mv = mv

        return self.mv

    def SingleView(self, model=sklearn.linear_model.LogisticRegression, model_params=None):
        concat_data = pd.concat(self.omics_data.values(), axis=1)
        self.sspa_scores_sv = self.sspa_method(concat_data, self.pathway_source, self.min_coverage)

        if model_params:
            sv = model(**model_params)
        else:
            sv = model()
        sv.fit(X=self.sspa_scores_sv, y=self.labels)
        sv.sspa_scores = self.sspa_scores_sv
        sv.name = 'SingleView'
        sv.coverage = self.coverage
        self.sv = sv

        return self.sv
    
    # cross-validation approaches

    def MultiViewCV(self):
        pass

def VIP_multiBlock(x_weights, x_superscores, x_loadings, y_loadings):
    # stack the weights from all blocks 
    weights = np.vstack(x_weights)
    # normalise the weights
    weights_norm = weights / np.sqrt(np.sum(weights**2, axis=0))
    # calculate product of sum of squares of superscores and y loadings
    sumsquares = np.sum(x_superscores**2, axis=0) * np.sum(y_loadings**2, axis=0)
    # p = number of variables - stack the loadings from all blocks
    p = np.vstack(x_loadings).shape[0]
    
    # VIP is a weighted sum of squares of PLS weights 
    vip_scores = np.sqrt(p * np.sum(sumsquares*(weights_norm**2), axis=1) / np.sum(sumsquares))
    return vip_scores

metab = pd.read_csv('data/metabolomics_example.csv', index_col=0)
prot = pd.read_csv('data/proteomics_example.csv', index_col=0)

# make possible to download MO paths from reactome
# mo_paths = sspa.process_reactome(
#     organism='Homo sapiens',
#     download_latest=True,
#     omics_type='multiomics',
#     filepath='data/')

# load pre-loaded pathways 
mo_paths = sspa.process_gmt(infile='data/Reactome_Homo_sapiens_pathways_multiomics_R85.gmt')

pi_model = PathIntegrate({'Metabolomics': metab, 'Proteomics':prot.iloc[:, :-1]}, metadata=prot['Group'], pathway_source=mo_paths, sspa_scoring='svd')

covid_multi_view = pi_model.MultiView(ncomp=5)

# launch the pathwy network explorer on a local server
launch_network_app(covid_multi_view, mo_paths)

# print(covid_multi_view.A_corrected_)
# print(covid_multi_view.vip)

# plot_functs.plot_block_importance(covid_multi_view)

# covid_single_view = pi_model.SingleView(model_params={'random_state':0})
# print(covid_single_view.intercept_)