from setuptools import setup
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()


setup(
    name='PathIntegrate',
    version='0.0.3',
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
        'pandas',
        'scikit-learn',
        'cmcrameri',
        'numpy',
        'setuptools',
        'scipy',
        'statsmodels',
        'sspa',
        'matplotlib',
        'seaborn',
        'dash_cytoscape',
        'plotly',
        'networkx',
        'dash_bootstrap_components',
        'dash'
    ],
)
