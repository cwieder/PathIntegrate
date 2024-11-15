# PathIntegrate
PathIntegrate Python package for pathway-based multi-omics data integration

![PyPI - Version](https://img.shields.io/pypi/v/PathIntegrate)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14060723.svg)](https://doi.org/10.5281/zenodo.14060723)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![PyPI - Downloads](https://img.shields.io/pypi/dm/pathintegrate)


![PathIntegrate graphical abstract](ModellingFrameworks_white.png "PathIntegrate graphical abstract")

#### Abstract
>As terabytes of multi-omics data are being generated, there is an ever-increasing need for methods facilitating the integration and interpretation of such data. Current multi-omics integration methods typically output lists, clusters, or subnetworks of molecules related to an outcome. Even with expert domain knowledge, discerning the biological processes involved is a time-consuming activity. Here we propose PathIntegrate, a method for integrating multi-omics datasets based on pathways, designed to exploit knowledge of biological systems and thus provide interpretable models for such studies. PathIntegrate employs single-sample pathway analysis to transform multi-omics datasets from the molecular to the pathway-level, and applies a predictive single-view or multi-view model to integrate the data. Model outputs include multi-omics pathways ranked by their contribution to the outcome prediction, the contribution of each omics layer, and the importance of each molecule in a pathway. 

## Features
- Pathway-based multi-omics data integration using PathIntegrate Multi-View and Single-View models
    - Multi-View model: Integrates multiple omics datasets using a shared pathway-based latent space
    - Single-View model: Integrates multi-omics data into one set of multi-omics pathway scores and applies an SKlearn-compatible predictive model
    - Pathway importance
    - Sample prediction
    - NEW unsupervised SingleView models (dimensionality reduction and clustering in the pathway space)
- SKlearn-like API for easy integration into existing pipelines
- Support for multiple pathway databases, including KEGG, Reactome, PathBank, and custom GMT files 
- Support for multiple pathway scoring methods available via the [sspa](https://github.com/cwieder/py-ssPA) package
- Cytoscape Network Viewer app for visualizing pathway-based multi-omics data integration results

![PathIntegrate Cytoscape app](App_network_view.JPG "Network viewer")

## System requirements
### Hardware
At least 8BG RAM recommended. PathIntegrate models can run on a Google Colab notebook server (see walkthrough tutorial below with example data).

### Software
PathIntegrate has been tested on MacOs, Windows 10 and Linux. Python 3.10 or higher is required. Python dependencies are listed in the requirements.txt file.

## Installation
```bash
pip install PathIntegrate
```

## Tutorials and documentation
Please see our Quickstart guide on [Google Colab](https://colab.research.google.com/drive/1MmGJp8I4JaIgFGNihYjaa4KJqYTGzdtv?usp=sharing)

Full documentation and function reference for PathIntegrate can be found via our [ReadTheDocs page](https://cwieder.github.io/PathIntegrate/)

## Citing PathIntegrate
If you use PathIntegrate in your research, please consider citing our paper:
```bibtex
@article{Wieder2024,
   author = {Cecilia Wieder and Juliette Cooke and Clement Frainay and Nathalie Poupin and Russell Bowler and Fabien Jourdan and Katerina J. Kechris and Rachel P.J. Lai and Timothy Ebbels},
   doi = {10.1371/JOURNAL.PCBI.1011814},
   issue = {3},
   journal = {PLOS Computational Biology},
   month = {3},
   pages = {e1011814},
   pmid = {38527092},
   publisher = {Public Library of Science},
   title = {PathIntegrate: Multivariate modelling approaches for pathway-based multi-omics data integration},
   volume = {20},
   url = {https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011814},
   year = {2024},
}
```

## PathIntegrate applications
Check out the following papers to see how PathIntegrate has been used in research:
- [Time-resolved, integrated multi-omic analysis reveals central role of amino acid pathways for defense responses in Arabidopsis thaliana](https://www.biorxiv.org/content/10.1101/2024.08.27.609849v1.full)

## License
GNU GPL v3

## Contributors
- Jude Popham @[judepops](https://github.com/judepops)