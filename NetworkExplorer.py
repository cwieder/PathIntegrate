# PathIntegrate network explorer app
# Takes a PathIntegrate model as input

# import packages
import pandas as pd
import numpy as np  
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

# Load extra layouts
cyto.load_extra_layouts()

# Set stylesheet
app = Dash(__name__, external_stylesheets=[dbc.themes.PULSE])

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
root_cmap = dict(zip(set(hierarchy_hsa_all['Root']), sns.color_palette("hls", len(set(hierarchy_hsa_all['Root']))).as_hex()))
cy_mo = nx.readwrite.json_graph.cytoscape_data(G)


def create_network(pi_model, pathway_source, hierarchy_source='preloaded'):
    pathways_accessible = pi_model.

    # Add model attributes to network
    name_dict = dict(zip(pathway_source.index, pathway_source['Pathway_name']))
    G.add_nodes_from([(node, {'Name': attr, 'label': attr}) for (node, attr) in name_dict.items()])
    G.add_nodes_from([(node, {'Root': attr, 
                              'RootCol': root_cmap[attr], 
                              'color': root_cmap[attr], 
                              'RootName': name_dict[attr]}) for (node, attr) in dict(zip(hierarchy_hsa_all[1], hierarchy_hsa_all['Root'])).items()])


    # add beta as node colour
    betas_cmap = dict(zip(all_pathways, sns.color_palette("RdBu", len(all_pathways)).as_hex()))
    G.add_nodes_from([(node, {'BetaColour': attr}) for (node, attr) in betas_cmap.items()])

    # add vip as node colour
    vip_cmap = dict(zip(all_pathways, get_hex_colors(vip_df['VIP_scaled'].tolist(), 'Blues')))
    G.add_nodes_from([(node, {'VIPColour': attr}) for (node, attr) in vip_cmap.items()])
    # start local server
    app.run_server(debug=True, use_reloader=False)




# Network layout configuration
SIDEBAR_STYLE = {
    "position": "fixed",
    "top": 0,
    "left": 0,
    "bottom": 0,
    "width": "16rem",
    "padding": "2rem 1rem",
    "background-color": "#BBDEFB",
    "padding": "6rem 1rem 2rem",
}
CONTENT_STYLE = {
    "margin-left": "18rem",
    "margin-right": "18rem",
    "padding": "2rem 1rem",
    # "border":"2px black solid"
    "background-color": "#FFFFFF",
}

sidebar = html.Div(
    [
    html.P("Data input"),
    html.Hr(),
    dcc.Upload(
        id='upload-data',
        children=html.Div([
            'Select files'
        ]),
        style={
            'width': '90%',
            'height': '60px',
            'lineHeight': '60px',
            'borderWidth': '1px',
            'borderStyle': 'dashed',
            'borderRadius': '5px',
            'textAlign': 'center',
            'margin': '5px'
        },
        # Allow multiple files to be uploaded
        multiple=True
    ),
    html.Div(id='output-data-upload'),

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
        for name in ['random', 'cose', 'circle', 'grid', 'concentric', 'cola', 'cise']
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
        for name in ['hierarchy', 'beta', 'VIP', 'P-value', 'pathway score']
    ]),
    html.P(
        "Show omics"
    ),
    dcc.Dropdown(
    id='dropdown-update-omics',
    value='metabolomics',
    clearable=False,
    options=[
        {'label': name.capitalize(), 'value': name}
        for name in ['metabolomics', 'proteomics', 'multi-omics']
    ]),
    html.P(
        "Export image"
    ),
    # dcc.Dropdown(
    # id='dropdown-export-img',
    # value='PNG',
    # clearable=False,
    # options=[
    #     {'label': name.upper(), 'value': name}
    #     for name in ['png', 'svg']
    # ]),
    html.Br(),
    dbc.Button("SVG", color="primary", id='btn-get-svg'),
    dbc.Button("PNG", color="primary", id='btn-get-png'),
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
    dbc.Input(id="input-pathway", placeholder="Select a pathway...", type="text"),
    dcc.Graph(id="bar-plot"),
    html.Br(),

    ],
    style={
    "position": "fixed",
    "top": 0,
    "right": 0,
    "bottom": 0,
    "width": "20rem",
    "padding": "2rem 1rem",
    "background-color": "#BBDEFB",
    "padding": "6rem 1rem 2rem",
},
)


navbar = dbc.NavbarSimple(
    children=[
        dbc.NavItem(dbc.NavLink("Page 1", href="#")),
        dbc.DropdownMenu(
            children=[
                dbc.DropdownMenuItem("More pages", header=True),
                dbc.DropdownMenuItem("Page 2", href="#"),
                dbc.DropdownMenuItem("Page 3", href="#"),
            ],
            nav=True,
            in_navbar=True,
            label="More",
        ),
    ],
    brand="Multi-omics pathway network explorer",
    brand_href="#",
    color="primary",
    dark=True,
    fixed='top'
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
            'width': '10',
            'height':'10',
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
content = html.Div(
    cyto.Cytoscape(
        id='mo_graph',
        layout={'name': 'random'},
        style={'width': '100%', 'height': '800px'},
        elements=
            cy_mo['elements']['nodes'] + cy_mo['elements']['edges'],
        stylesheet=default_stylesheet
    ), style=CONTENT_STYLE
)

app.layout = html.Div([dcc.Location(id="url"), navbar, sidebar, content, sidebar2])

# Callback for updating the graph layout

# callbacks
@app.callback(Output('mo_graph', 'layout'), Input('dropdown-update-layout', 'value'))
def update_layout(layout):
    return {
        'name': layout,
        'animate': False
    }



@app.callback(Output('mo_graph', 'stylesheet'),
              Input('dropdown-update-node-color', 'value'), Input('dropdown-update-omics', 'value'))
def update_stylesheet(node_value, omics):
    node_value_color = {'hierarchy': 'data(color)', 'beta': 'data(BetaColour)', 'VIP': 'data(VIPColour)', 'P-value': 'data(PColour)'}
    omics_selector = {'metabolomics': 'data(Coverage_m)', 'proteomics': 'data(Coverage_p)', 'multi-omics': 'data(Coverage_mo)'}

    if node_value != 'pathway score':
        new_styles = [
        {
            'selector': 'node',
            'style': {
                'background-color': node_value_color[node_value],
                'border-width': '1px',
                'border-color': 'black',
                'background-fit': 'cover',
                'height': omics_selector[omics],
                'width': omics_selector[omics]
                # 'background-image': bg_img[node_value]
            }
        },
        {
            'selector': 'edge',
            'style': {
                'line-color': '#A3C4BC'
            }
        }]

    else:
        new_styles = [
        {
            'selector': 'node',
            'style': {
                'border-width': '1px',
                'border-color': 'black',
                'background-fit': 'cover',
                'background-image': 'data(NodeSVG)',
                'height': omics_selector[omics]**10,
                'width': omics_selector[omics]**10
            }
        },
        {
            'selector': 'edge',
            'style': {
                'line-color': '#A3C4BC'
            }
        }]

    return default_stylesheet + new_styles

@app.callback(Output('mo_graph', 'elements'), Input('dropdown-update-omics', 'value'))
def update_layout(omics):
    omics_selector = {'multi-omics': cy_mo}
        # omics_selector = {'metabolomics': cy, 'proteomics': cy_p, 'multi-omics': cy_mo}
    return omics_selector[omics]['elements']['nodes'] + omics_selector[omics]['elements']['edges']

# @app.callback(Output('mo_graph', 'stylesheet'), Input('dropdown-update-omics', 'value'))
# def update_size(omics):
#     omics_selector = {'metabolomics': data['Coverage_m'], 'proteomics': data['Coverage_P']}

#     new_styles = [
#         {
#             'selector': 'node',
#             'style': {
#                 'background-color': node_value_color[node_value],
#                 'border-width': '1px',
#                 'border-color': 'black',
#                 'background-fit': 'cover',
#                 'height': omics_selector[omics],
#                 'width': omics_selector[omics]
#                 # 'background-image': bg_img[node_value]
#             }
#         }]
#     return default_stylesheet + new_styles

# show hover text
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
        cvrg1 = None
        cvrg2 = None
        try:
            cvrg1 = data['Coverage_m']
        except KeyError:
            cvrg1 = 'Not present in M' 

        try:
            cvrg2 = data['Coverage_p']
        except KeyError:
            cvrg2 = 'Not present in P' 


        return str(cvrg1) + " " + str(cvrg2)

# Get user input files
def parse_contents(contents, filename, date):
    content_type, content_string = contents.split(',')

    decoded = base64.b64decode(content_string)
    try:
        if 'csv' in filename:
            # Assume that the user uploaded a CSV file
            df = pd.read_csv(
                io.StringIO(decoded.decode('utf-8')))
        elif 'xls' in filename:
            # Assume that the user uploaded an excel file
            df = pd.read_excel(io.BytesIO(decoded))
        return df
    except Exception as e:
        print(e)
        return html.Div([
            'There was an error processing this file.'
        ])
    

@app.callback(Output('output-data-upload', 'children'),
              Input('upload-data', 'contents'),
              State('upload-data', 'filename'),
              State('upload-data', 'last_modified'))
def update_output(list_of_contents, list_of_names, list_of_dates):
    if list_of_contents is not None:
        children = list_of_names
            # parse_contents(c, n, d) for c, n, d in
            # zip(list_of_contents, list_of_names, list_of_dates)]
        return children

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
