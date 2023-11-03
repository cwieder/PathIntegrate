from setuptools import setup
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='PathIntegrate',
    version='0.0.1',
    packages=['pathintegrate'],
    package_dir={'':'src'},
    url='https://github.com/cwieder/PathIntegrate',
    license='GNU 3.0',
    author='Cecilia Wieder',
    author_email='cw2019@ic.ac.uk',
    description='PathIntegrate: multivariate modelling approaches for pathway-based muti-omics integration',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={'pathintegrate': ['data/*']},
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
    ],
    install_requires=[
        'cmcrameri==1.3',
        'dash==2.14.1',
        'dash_bootstrap_components==1.5.0',
        'dash_cytoscape==0.3.0',
        'datauri==1.0.0',
        'matplotlib==3.8.1',
        'mbpls==1.0.4',
        'networkx==3.2.1',
        'numpy==1.26.1',
        'pandas==2.1.2',
        'plotly==5.18.0',
        'scikit_learn==1.3.2',
        'scipy==1.11.3',
        'seaborn==0.13.0',
        'setuptools==68.0.0',
        'sspa>=1.0.1',
        'statsmodels==0.14.0',
        'svgwrite==1.4.3'
    ],
)
