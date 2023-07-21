# import packages
import pandas as pd
import numpy as np  
import dash
from dash import Dash, html, dcc, Input, Output, State, ctx
import plotly.express as px
import dash_cytoscape as cyto
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import networkx as nx
import dash_bootstrap_components as dbc
import matplotlib.cm as cm
import matplotlib as matplotlib
import svgwrite
from datauri import DataURI
from dash_bootstrap_components._components.Container import Container
from pathlib import Path
dash.register_page(__name__, path='/details')


layout = html.Div(children=[
    html.H1(children='Pathway details view'),

    html.Div(children='''
        This is our Home page content.
    '''),

])