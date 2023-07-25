# PathIntegrate network explorer app
# Takes a PathIntegrate model as input

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
from dash.exceptions import PreventUpdate
from plotly.subplots import make_subplots
import plotly.graph_objects as go


# Save downloads to default Downloads folder
downloads_path = str(Path.home() / "Downloads")

# Load extra layouts
cyto.load_extra_layouts()

# Set stylesheet
app = Dash(__name__, external_stylesheets=[dbc.themes.YETI], use_pages=False)

# generate node colours based on model attributes
def get_hex_colors(values, cmap):
    norm = plt.Normalize(min(values), max(values))
    cmap = plt.cm.get_cmap(cmap)  # Change the colormap here as desired

    hex_colors = [mcolors.rgb2hex(cmap(norm(value))) for value in values]
    return hex_colors

# find parent pathway for each pathway
def find_root(G,child):
    parent = list(G.predecessors(child))
    if len(parent) == 0:
        return child
    else:  
        return find_root(G, parent[0])
    
# load the pathway database file from the data folder
hierarchy = pd.read_csv('data/ReactomePathwaysRelation.txt', sep='\t', header=None)
hierarchy_hsa = hierarchy[hierarchy[0].str.contains('HSA')]
hierarchy_hsa_parents = np.setdiff1d(hierarchy_hsa[0], hierarchy_hsa[1])
hierarchy_hsa_all = pd.concat([hierarchy_hsa, pd.DataFrame([hierarchy_hsa_parents, hierarchy_hsa_parents], index=[0, 1]).T])

# the default graph is the pathway hierarchy coloured by root pathway membership as defined by Reactome
G = nx.from_pandas_edgelist(hierarchy_hsa, source=0, target=1, create_using=nx.DiGraph())
hierarchy_hsa_all['Root'] = [find_root(G, i) for i in hierarchy_hsa_all[1]]
root_cmap = dict(zip(set(hierarchy_hsa_all['Root']), sns.color_palette("husl", len(set(hierarchy_hsa_all['Root']))).as_hex()))
cy_mo = nx.readwrite.json_graph.cytoscape_data(G)



# Network layout configuration
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "1rem",
    "padding-top": "5rem",
    "background-color": "#BBDEFB",
    # "padding": "6rem 1rem 2rem",
}
CONTENT_STYLE = {
    "margin-left": "16rem",
    "margin-right": "16rem",
    # "padding": "1rem 1rem",
    "padding-top": "4rem",
    "padding-bottom": "1rem",
    "padding-left": "1rem",
    "padding-right": "1rem",
    # "border":"2px black solid"
    "background-color": "#FFFFFF",
}

sidebar = html.Div(
    [

    html.P("Network options"),
    html.Hr(),

    html.P(
        "Layout option"
    ),
    dcc.Dropdown(
    id='dropdown-update-layout',
    value='random',
    clearable=False,
    options=[
        {'label': name.capitalize(), 'value': name}
        for name in ['random', 'cose', 'circle', 'grid', 'concentric', 'cola']
    ]),
    html.P(
        "Colour nodes by"
    ),
    dcc.Dropdown(
    id='dropdown-update-node-color',
    value='hierarchy',
    clearable=False,
    options=[
        {'label': name.capitalize(), 'value': name}
        for name in ['hierarchy', 'beta', 'VIP']
    ]),
    # html.P(
    #     "Show omics"
    # ),
    # dcc.Dropdown(
    # id='dropdown-update-omics',
    # value='multi-omics',
    # clearable=False,
    # options=[
    #     {'label': name.capitalize(), 'value': name}
    #     for name in ['metabolomics', 'proteomics', 'multi-omics']
    # ]),
    html.P(
        "Export network"
    ),

    dbc.Button("SVG", color="primary", id='btn-get-svg'),
    dbc.Button("PNG", color="primary", id='btn-get-png'),
    dbc.Button("Network", color="primary", id='btn-get-gml'),
    dcc.Download(id="download-network")
    ],
    style=SIDEBAR_STYLE,
)


sidebar2 = html.Div(
    [html.P("Node information"),
    html.Hr(),
    dbc.ListGroup(
    [
        dbc.ListGroupItem(
            html.Div(
                [html.P("Pathway name"), html.P(id='cytoscape-mouseoverNodeData-output-name')
                ])),
        dbc.ListGroupItem(html.Div(
                [html.P("Parent pathway"), html.P(id='cytoscape-mouseoverNodeData-output-root')
                ])),
        dbc.ListGroupItem(html.Div(
                [html.P("Coverage"), html.P(id='cytoscape-mouseoverNodeData-output-coverage')
                ])),
    ]),

    html.Br(),

    ],
    style={
    "position": "fixed",
    "top": 0,
    "right": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "1rem",
    "padding-top": "5rem",
    "background-color": "#BBDEFB",
},
)



navbar = dbc.NavbarSimple(
    # children=[
    #     dbc.NavItem(dbc.NavLink("Home", href="/")),
    #     dbc.DropdownMenu(
    #         children=[
    #             dbc.DropdownMenuItem("Detail view", href="/details"),
    #             dbc.DropdownMenuItem("About", href="/about"),

    #         ],
    #         nav=True,
    #         in_navbar=True,
    #         label="More",
    #     ),
    # ],
    brand="PathIntegrate multi-omics pathway network explorer",
    brand_href="#",
    color="primary",
    dark=True,
    fixed='top',
    style={
    # "position": "fixed",
    # "top": 0,
    # "right": 0,
    # "bottom": 0,
    # "padding-bottom": "2rem",
    # "width": "20rem",
    # "padding": "2rem 1rem",
    # "background-color": "#BBDEFB",
    # "padding": "6rem 1rem 2rem",
    },
)

# default stylesheet
default_stylesheet = [
    {
        'selector': 'node',
        'style': {
            'background-color': 'data(color)',
            'shape': 'ellipse',
            'label': 'data(label)',
            'text-wrap': 'wrap',
            'text-background-color': 'yellow',
            'text-max-width': '120px',
            'width': 'data(MO_coverage)',
            'height':'data(MO_coverage)',
            'text-justification': 'auto',
            'font-family': ['Verdana', 'Roboto', 'Arial'],
            'font-size': '10px'
        }
    },
    {
        'selector': 'edge',
        'style': {
            'line-color': '#A3C4BC'
        }
    }
]



# Callback for updating the graph layout
@app.callback(Output('mo_graph', 'layout'), Input('dropdown-update-layout', 'value'))
def update_layout(layout):
    return {
        'name': layout,
        'animate': True
    }

# callback for updating node colour
@app.callback(Output('mo_graph', 'stylesheet'),
              Input('dropdown-update-node-color', 'value'))
def update_stylesheet(node_value):
    node_value_color = {'hierarchy': 'data(color)', 'beta': 'data(BetaColour)', 'VIP': 'data(VIPColour)'}

    new_styles = [
        {
            'selector': 'node',
            'style': {
                'background-color': node_value_color[node_value],
                'border-width': '1px',
                'border-color': 'black',
                'background-fit': 'cover',
                'height': 'data(MO_coverage)',
                'width':'data(MO_coverage)'
                # 'background-image': bg_img[node_value]
            }
        },
        {
            'selector': 'edge',
            'style': {
                'line-color': '#A3C4BC'
            }
    }]


    return default_stylesheet + new_styles


@app.callback(Output('cytoscape-mouseoverNodeData-output-name', 'children'),
              Input('mo_graph', 'mouseoverNodeData'))
def displayTapNodeData(data):
    if data:
        return data['label']

@app.callback(Output('cytoscape-mouseoverNodeData-output-root', 'children'),
              Input('mo_graph', 'mouseoverNodeData'))
def displayTapNodeData(data):
    if data:
        return data['Root']

@app.callback(Output('cytoscape-mouseoverNodeData-output-coverage', 'children'),
              Input('mo_graph', 'mouseoverNodeData'))
def displayTapNodeData(data):
    if data:
        return data['MO_coverage']


# Download image
@app.callback(
    Output("mo_graph", "generateImage"),
    [
        Input("btn-get-svg", "n_clicks"),
        Input("btn-get-png", "n_clicks")
    ])
def get_image(get_clicks_svg, get_clicks_png):

    # File type to output of 'svg, 'png', 'jpg', or 'jpeg' (alias of 'jpg')

    # 'store': Stores the image data in 'imageData' !only jpg/png are supported
    # 'download'`: Downloads the image as a file with all data handling
    # 'both'`: Stores image data and downloads image as file.
    ftype='png'
    action = 'store'

    if ctx.triggered:
        if ctx.triggered_id != "tabs":
            action = "download"
            ftype = ctx.triggered_id.split("-")[-1]
            print(ftype)

    return {
        'type': ftype,
        'action': action
        }

@app.callback(
    Output("download-network", "data"),
    Input("btn-get-gml", "n_clicks"),
    prevent_initial_call=True,
)
def download_network(n_clicks):
    download_path = downloads_path+"/PathIntegrate_network.gml"
    nx.write_gml(MO_graph, download_path)
    print("Network .gml file saved to: "+download_path)

    return nx.write_gml(MO_graph, download_path)

    
@app.callback([Output("fig_molecular", "figure"),
                Output("pathway_selected", "children")],
                  Input("dropdown", "value"))
def update_bar_chart(pathway):
    if modelname == 'MultiView':
        pathways_dfs = []
        for k, v in molecule_importances.items():
            try:
                pdf = v[pathway]
                pdf = pdf.add_suffix(k)
                pathways_dfs.append(pdf)
            except KeyError:    
                pass

        pathway_df_molec = pd.concat(pathways_dfs, axis=1)
        fig = make_subplots(rows=1, cols=len(pathways_dfs), shared_xaxes='rows')

        for i in range(0, len(pathways_dfs)):
            fig.add_trace(go.Bar(x=pathway_df_molec.index, y=pathway_df_molec.iloc[:,i], name=pathway_df_molec.columns[i]), row=1, col=i+1)
 
        return fig, name_dict[pathway]
    else:
        pathway_df_molec = molecule_importances[pathway]
        fig = go.Figure()
        fig.add_trace(go.Bar(x=pathway_df_molec.index.tolist(), y=pathway_df_molec['loadings']))
        return fig, name_dict[pathway]
# molecule level vis
# @app.callback(
#     Output("bar-plot", "figure"), 
#     Input("input-pathway", "value"))
# def update_bar_chart(data):
#     if data:
#         mean_vals = metab.loc[:, metab.columns.isin(mo_paths_dict[data] + ['Group'])].groupby('Group').mean()
#         mean_vals_long = mean_vals.melt(ignore_index=False).reset_index()
#         fig = px.bar(mean_vals_long, x="variable", y="value", 
#                     color="Group", barmode="group")
#         return fig


# start local server
def launch_network_app(pi_model, pathway_source, hierarchy_source='preloaded'):
    global pathways_accessible

    # Add model attributes to network
    global name_dict
    name_dict = dict(zip(pathway_source.index, pathway_source['Pathway_name']))
    G.add_nodes_from([(node, {'Name': attr, 'label': attr}) for (node, attr) in name_dict.items()])
    G.add_nodes_from([(node, {'Root': attr, 
                              'RootCol': root_cmap[attr], 
                              'color': root_cmap[attr], 
                              'RootName': name_dict[attr]}) for (node, attr) in dict(zip(hierarchy_hsa_all[1], hierarchy_hsa_all['Root'])).items()])
    G.add_nodes_from([(node, {'MO_coverage': attr}) for (node, attr) in pi_model.coverage.items()])

    global modelname
    modelname = pi_model.name

    if pi_model.name == 'MultiView':
        pathways_accessible = list(set(sum([i.columns.tolist() for i in pi_model.sspa_scores.values()], [])))
        # add beta as node colour
        betas_cmap = dict(zip(pathways_accessible, get_hex_colors(pi_model.beta, 'RdBu')))
        G.add_nodes_from([(node, {'BetaColour': attr}) for (node, attr) in betas_cmap.items()])

        # add vip as node colour
        vip_cmap = dict(zip(pathways_accessible, get_hex_colors(pi_model.vip['VIP_scaled'].tolist(), 'Blues')))
        G.add_nodes_from([(node, {'VIPColour': attr}) for (node, attr) in vip_cmap.items()])

    if pi_model.name == 'SingleView':
        pathways_accessible = pi_model.sspa_scores.columns.tolist()
        # add beta as node colour
        # betas_cmap = dict(zip(pathways_accessible, get_hex_colors(pi_model.beta, 'RdBu')))
        # G.add_nodes_from([(node, {'BetaColour': attr}) for (node, attr) in betas_cmap.items()])

        # # add vip as node colour
        # vip_cmap = dict(zip(pathways_accessible, get_hex_colors(pi_model.vip['VIP_scaled'].tolist(), 'Blues')))
        # G.add_nodes_from([(node, {'VIPColour': attr}) for (node, attr) in vip_cmap.items()])

    # add molecular importances for plotting
    global molecule_importances
    molecule_importances = pi_model.molecular_importances

    
    # only show nodes with sufficient coverage
    global MO_graph
    MO_graph = G.subgraph(pathways_accessible)
    cy_mo = nx.readwrite.json_graph.cytoscape_data(MO_graph)
    # network generation
    content = html.Div(
        cyto.Cytoscape(
            id='mo_graph',
            layout={'name': 'random'},
            style={'width': '100%', 'height': '800px'},
            elements=cy_mo['elements']['nodes'] + cy_mo['elements']['edges'],
            stylesheet=default_stylesheet
        ),
    # style=CONTENT_STYLE
    )

    app.layout = html.Div([
        navbar,
        sidebar,
        sidebar2,
        dcc.Tabs([
            dcc.Tab(label='Network', children=[
                content,
            ]),
            dcc.Tab(label='Molecular importance', children=[
                html.Div([
                    html.Br(),
                    html.H4(
                    "Select a pathway from the dropdown menu to view molecular importance:"),
                    dcc.Dropdown(pathways_accessible,
                                 pathways_accessible[0],
                                 placeholder="Select a pathway",
                    id="dropdown"),
                    html.Br(),
                    html.P(id="pathway_selected")
                ]),
                dcc.Graph(id="fig_molecular"),
                
            ]),
        ])
        ],
        style=CONTENT_STYLE

    )

    # app.layout = dbc.Container([
    #                         html.Div([ dcc.Location(id="url",refresh=False), 
    #                         navbar, 
    #                         sidebar,
    #                         content,
    #                         sidebar2,
    #                         ]),],fluid=True)
    # app.layout = html.Div([dcc.Location(id="url"), navbar, sidebar, content, sidebar2])
    app.run_server(debug=True, use_reloader=True)


 
   