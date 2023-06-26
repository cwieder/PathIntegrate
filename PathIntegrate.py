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
        mv.fit(sspa_scores, self.labels)

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


metab = pd.read_csv('metabolomics_example.csv', index_col=0)
prot = pd.read_csv('proteomics_example.csv', index_col=0)

# make possible to download MO paths from reactome
mo_paths = pd.read_csv("../Pathway_databases/Reactome_multi_omics_ChEBI_Uniprot.csv", index_col=0)

pi_model  = PathIntegrate({'Metabolomics': metab, 'Proteomics':prot.iloc[:, :-1]}, metadata=prot['Group'], pathway_source=mo_paths, sspa_scoring='zscore')

covid_multi_view = pi_model.MultiView(ncomp=5)
print(covid_multi_view.A_corrected_)

covid_single_view = pi_model.SingleView(model_params={'random_state':0})
print(covid_single_view.intercept_)