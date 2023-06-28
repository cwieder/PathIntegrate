import pandas as pd
import numpy as np
import sspa
import sklearn
from mbpls.mbpls import MBPLS

class PathIntegrate:

    def __init__(self, omics_data:dict, metadata, pathway_source, sspa_scoring='svd'):
        self.omics_data = omics_data
        self.metadata = metadata
        self.pathway_source = pathway_source
        self.sspa_scoring = sspa_scoring
        
        sspa_methods = {'svd': sspa.sspa_svd, 'ssGSEA': sspa.sspa_ssGSEA, 'kpca': sspa.sspa_kpca, 'ssClustPA': sspa.sspa_ssClustPA, 'zscore': sspa.sspa_zscore}
        self.sspa_method = sspa_methods[self.sspa_scoring]

        self.labels = pd.factorize(self.metadata)[0]

    def MultiView(self, ncomp=2):
        sspa_scores = [self.sspa_method(i, self.pathway_source) for i in self.omics_data.values()]
        mv = MBPLS(n_components=ncomp)
        mv.fit([i.copy(deep=True) for i in sspa_scores], self.labels)
        vip_scores = VIP_multiBlock(mv.W_, mv.Ts_, mv.P_, mv.V_)
        vip_df = pd.DataFrame(vip_scores, index=sum([i.columns.tolist() for i in sspa_scores], []))
        vip_df['Name'] = vip_df.index.map(dict(zip(self.pathway_source.index, self.pathway_source['Pathway_name'])))
        mv.vip = vip_df
        return mv

    def SingleView(self, model=sklearn.linear_model.LogisticRegression, model_params=None):
        concat_data = pd.concat(self.omics_data.values(), axis=1)
        sspa_scores = self.sspa_method(concat_data, self.pathway_source)

        if model_params:
            sv = model(**model_params)
        else:
            sv = model()
        sv.fit(X=sspa_scores, y=self.labels)

        return sv
    
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

metab = pd.read_csv('metabolomics_example.csv', index_col=0)
prot = pd.read_csv('proteomics_example.csv', index_col=0)

# make possible to download MO paths from reactome
mo_paths = pd.read_csv("../Pathway_databases/Reactome_multi_omics_ChEBI_Uniprot.csv", index_col=0)

pi_model  = PathIntegrate({'Metabolomics': metab, 'Proteomics':prot.iloc[:, :-1]}, metadata=prot['Group'], pathway_source=mo_paths, sspa_scoring='zscore')

covid_multi_view = pi_model.MultiView(ncomp=5)
print(covid_multi_view.A_corrected_)
print(covid_multi_view.vip)

covid_single_view = pi_model.SingleView(model_params={'random_state':0})
print(covid_single_view.intercept_)