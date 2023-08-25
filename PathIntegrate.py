import pandas as pd
import numpy as np
import sspa
import sklearn
from mbpls.mbpls import MBPLS
import plot_functs
from app import launch_network_app
from sklearn.preprocessing import StandardScaler

class PathIntegrate:

    def __init__(self, omics_data:dict, metadata, pathway_source, sspa_scoring='svd', min_coverage=3):
        self.omics_data = omics_data
        self.metadata = metadata
        self.pathway_source = pathway_source
        self.pathway_dict = sspa.utils.pathwaydf_to_dict(pathway_source)
        self.sspa_scoring = sspa_scoring
        self.min_coverage = min_coverage

        sspa_methods = {'svd': sspa.sspa_SVD, 'ssGSEA': sspa.sspa_ssGSEA, 'kpca': sspa.sspa_KPCA, 'ssClustPA': sspa.sspa_ssClustPA, 'zscore': sspa.sspa_zscore}
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
        print('Generating pathway scores...')
        sspa_scores = [self.sspa_method(self.pathway_source, self.min_coverage).fit_transform(i) for i in self.omics_data.values()]
        # sspa_scores = [self.sspa_method(i, self.pathway_source, self.min_coverage, return_molecular_importance=True) for i in self.omics_data.values()]

        self.sspa_scores_mv = dict(zip(self.omics_data.keys(), sspa_scores))
        print('Fitting MultiView model')
        mv = MBPLS(n_components=ncomp)
        mv.fit([i.copy(deep=True) for i in self.sspa_scores_mv.values()], self.labels)

        # compute VIP and scale VIP across omics
        vip_scores = VIP_multiBlock(mv.W_, mv.Ts_, mv.P_, mv.V_)
        vip_df = pd.DataFrame(vip_scores, index=sum([i.columns.tolist() for i in self.sspa_scores_mv.values()], []))
        vip_df['Name'] = vip_df.index.map(dict(zip(self.pathway_source.index, self.pathway_source['Pathway_name'])))
        vip_df['Source'] = sum([[k] * v.shape[1] for k, v in self.sspa_scores_mv.items()], [])
        vip_df['VIP_scaled'] = vip_df.groupby('Source')[0].transform(lambda x: StandardScaler().fit_transform(x.values[:,np.newaxis]).ravel())

        mv.name = 'MultiView'

        # only some sspa methods can return the molecular importance
        if hasattr(sspa_scores[0], 'molecular_importance'):
            mv.molecular_importances = dict(zip(self.omics_data.keys(), [i.molecular_importance for i in sspa_scores]))
        mv.beta = mv.beta_.flatten()
        mv.vip = vip_df
        mv.omics_names = list(self.omics_data.keys())
        mv.sspa_scores = self.sspa_scores_mv
        mv.coverage = self.coverage
        self.mv = mv

        return self.mv

    def SingleView(self, model=sklearn.linear_model.LogisticRegression, model_params=None):
        """Fits a PathIntegrate SingleView model using an SKLearn-compatible predictive model.

        Args:
            model (object, optional): SKlearn prediction model class. Defaults to sklearn.linear_model.LogisticRegression.
            model_params (_type_, optional): Model-specific hyperparameters. Defaults to None.

        Returns:
            object: Fitted PathIntegrate SingleView model.
        """
        concat_data = pd.concat(self.omics_data.values(), axis=1)
        print('Generating pathway scores...')

        sspa_scores = self.sspa_method(self.pathway_source, self.min_coverage).fit_transform(concat_data)
        self.sspa_scores_sv = sspa_scores
       
        if model_params:
            sv = model(**model_params)
        else:
            sv = model()
        print('Fitting SingleView model')

        sv.fit(X=self.sspa_scores_sv, y=self.labels)
        sv.sspa_scores = self.sspa_scores_sv
        sv.name = 'SingleView'
        sv.coverage = self.coverage

        # only some sspa methods can return the molecular importance
        if hasattr(sspa_scores, 'molecular_importance'):
            sv.molecular_importance = sspa_scores.molecular_importance
        self.sv = sv

        return self.sv
    
    # cross-validation approaches
    def SingleViewCV(self):
        pass

    def MultiViewCV(self):
        # Set up sklearn pipeline
        pipe_mv = sklearn.pipeline.Pipeline([
            ('mbpls', MBPLS(n_components=2))
        ])
        pass

        # Set up cross-validation



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

# metab = pd.read_csv('data/metabolomics_example.csv', index_col=0)
# prot = pd.read_csv('data/proteomics_example.csv', index_col=0)

# # make possible to download MO paths from reactome
# # mo_paths = sspa.process_reactome(
# #     organism='Homo sapiens',
# #     download_latest=True,
# #     omics_type='multiomics',
# #     filepath='data/')

# # load pre-loaded pathways 
# mo_paths = sspa.process_gmt(infile='data/Reactome_Homo_sapiens_pathways_multiomics_R85.gmt')

# pi_model = PathIntegrate({'Metabolomics': metab, 'Proteomics':prot.iloc[:, :-1]}, metadata=prot['Group'], pathway_source=mo_paths, sspa_scoring='svd', min_coverage=2)

# covid_multi_view = pi_model.MultiView(ncomp=5)

# # launch the pathwy network explorer on a local server
# launch_network_app(covid_multi_view, mo_paths)

# print(covid_multi_view.A_corrected_)
# print(covid_multi_view.vip)

# plot_functs.plot_block_importance(covid_multi_view)

# covid_single_view = pi_model.SingleView(model_params={'random_state':0})
# launch_network_app(covid_single_view, mo_paths)
# print(covid_single_view.intercept_)