import pandas as pd
import numpy as np
import pathintegrate
import sspa


md = pd.read_csv('H:/Documents/pathway-integration/COPDgene/COPDgene_phonotype.txt', sep='\t')
prot = pd.read_csv('H:/Documents/pathway-integration/COPDgene/COPDgene_proteomics_UniProt.csv', index_col=0)
metab = pd.read_csv('H:/Documents/pathway-integration/COPDgene/COPDgene_metabolomics_CHEBI_mapped.csv', index_col=0)
trans = pd.read_csv('D:/COPDgene/Processed/COPDgene_transcriptomics_filt_Q1_scaled.csv', index_col=0)
metab['Group'] = metab.index.map(dict(zip(md['sid'], md["COPD"])))

metab = metab[metab['Group'].isin([0, 1])]
intersect_samples = set(metab.index.tolist()) & set(prot.index.tolist()) & set(trans.index.tolist())
prot = prot.loc[intersect_samples, :]
metab = metab.loc[intersect_samples, :]
trans = trans.loc[intersect_samples, :]

mo_paths_all = pd.read_csv('D:/Pathway_databases\Reactome_multi_omics_ChEBI_Uniprot_Ensembl.csv', index_col=0, dtype=object)

pi_model = pathintegrate.PathIntegrate(omics_data={'Metabolomics': metab, 'Proteomics':prot.iloc[:, :-1], 'Transcriptomics': trans.iloc[:, :-1]}, 
                                       metadata=metab['Group'], 
                                       pathway_source=mo_paths_all, 
                                       sspa_scoring='svd', 
                                       min_coverage=2)

copdgene_multi_view = pi_model.MultiView(ncomp=4)

# launch the pathwy network explorer on a local server
pathintegrate.launch_network_app(copdgene_multi_view, mo_paths_all)