"""
;===================================================================================================
; Title:   Visualization of the acetylome changes in Pseudomonas aeruginosa upon phage infection
; Authors: Stefaan Verwimp, Aditya Badola, Hannelore Longin, Ben De Maesschalck
;===================================================================================================


"""
import pandas as pd
import dash	
from dash.dependencies import Input, Output, State	 	# this line will give an error if there is a file called 'dash.py' in the project
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table
import dash_cytoscape as cyto


#############################################################################################################
##	Variables that are free to be changed 																   ##
#############################################################################################################
																											#
# Cut-off score for combined_score from string interactions													#
NODE_CUTOFF_SCORE = 0.7			# Cut-off value for interactions between proteins							#
																											#
# Colors																									#
neutral_color = '#6c6f74'		# grey																		#
positive_color = '#7bb526'		# green																		#
negative_color = '#f32c22'		# red																		#
selected_edge_color = '#3c6975'	# blue-ish																	#
																											#
##############################################################################################################
colordict = {'positive': positive_color, 'negative': negative_color, 'similar': neutral_color}


# Read in the data
nodeDf = pd.read_csv('Preprocessing/Output/nodeDf.tsv', delimiter='\t')
prot_annot = pd.read_csv('Preprocessing/String_man/string_protein_annotations.tsv',sep='\t')
kegg = pd.read_csv('Preprocessing/Output/pathways.tsv', sep='\t')
acetylation = pd.read_csv('Preprocessing/Output/ackegg.tsv', sep='\t')

# The number of acetylation sites is read as floats by python, so we change them back into integers
acetylation['numAcSites'] = acetylation['numAcSites'].apply(lambda x: int(x))


# Splitting the different pathyways from the kegg file to make it more readable in the "selected node details" card in the application.
# The resulting list will be used in the "displaySelectedNodeData" callback.
kegg['keggPathways'] = kegg['keggPathways'].apply(lambda x: x.split(' // '))

# Dictionaries are for retrieving protein (node) information in a O(1) manner, in order to minimize delays.
protein_annotation = dict(zip(prot_annot['identifier'], prot_annot['annotation']))
kegg_dict = dict(zip(kegg['keggID'], kegg['keggPathways']))
prot_anno = dict(zip(prot_annot['node'],prot_annot['annotation']))
acet_sites = dict(zip(acetylation['uniprotID'], acetylation['numAcSites']))
logfc_anno = dict(zip(acetylation['uniprotID'], acetylation['protLogFC']))


#Make a list of unique records
unique_nodict = prot_anno
for k, v in unique_nodict.copy().items():
	if v == "annotation not available":
		del unique_nodict[k]
unique_keys = list(unique_nodict.keys())


# make nodes and edges
nodes = set()
cy_edges = []
cy_nodes = []


# remove all entries from dataframe that do not meet requirements for the interaction score
nodeDf = nodeDf[nodeDf.combined_score >= NODE_CUTOFF_SCORE]


# logfc_anno string is a bit messy, so this function returns a short, more orderly string
# something like: "negative !NaN in peptide(s)" will be returned as just "negative"
def get_logFC_as_string(logfc_anno):
	logfc = str(logfc_anno)
	if 'positive' in logfc:
		return 'positive'
	if 'negative' in logfc:
		return 'negative'
	else:
		return 'similar'

for index, row in nodeDf.iterrows():
	# source node annotation info
	source = row['node1']
	source_stringid = row['node1_string_id']
	source_uniprot = row['node1_uniprot']
	source_kegg = row['node1_kegg']
	source_logFC = get_logFC_as_string(logfc_anno.get(source_uniprot))

	# target node annotation info
	target = row['node2']
	target_stringid = row['node2_string_id']
	target_uniprot = row['node2_uniprot']
	target_kegg = row['node2_kegg']
	target_logFC = get_logFC_as_string(logfc_anno.get(target_uniprot))

	# edge annotation info
	score = row['combined_score']
	interaction = row['interaction']

	# put information for current row of the dataframe into dictionaries
	cy_source = {'data': {'id': source, 'label': source, 'stringid':source_stringid, 'UniprotID': source_uniprot, 'KEGG_ID': source_kegg, 'logFC': source_logFC}}
	cy_target = {'data': {'id': target, 'label': target, 'stringid':target_stringid, 'UniprotID': target_uniprot, 'KEGG_ID': target_kegg, 'logFC': target_logFC}}
	cy_edge = {'data': {'id': source+target, 'source': source, 'target': target, 'score': score, 'interaction': interaction, 'source_logFC': source_logFC, 'target_logFC': target_logFC}}

	# Add nodes if not already in the set
	if source not in nodes:
		nodes.add(source)
		cy_nodes.append(cy_source)	  # data: { [here comes that data contained for the node as a dictionary] }

	if target not in nodes:
		nodes.add(target)
		cy_nodes.append(cy_target)
	
	# No need to filter on combine_score since entries with score below 0.7 were removed before this
	cy_edges.append(cy_edge)


# Store original dataframes in different variable, so that this can be called when removing filtering in tabular view
nodeDf_orig = nodeDf
acetylation_orig = acetylation

# Default variable declared for label
label = 'data(label)'


# Default style of nodes and edges (will be loaded first)
default_stylesheet  = [
	{
		'selector': 'node',						# for all nodes
		'style': {
			'opacity': 1,
			'label': label,
			'background-color': neutral_color,	# background-color = color of node
			'color': '#000000',					# label color = black
		}
	},
	{
		'selector': '[logFC = "positive"]',		# now, only for nodes with a positive log fold change	(overwrites stylesheet for all nodes)
		'style': {
			'opacity': 1,
			'label': label,
			'background-color': positive_color,
		}

	},
	{
		'selector': '[logFC = "negative"]',		# now, only for nodes with a negative log fold change 	(overwrites stylesheet for all nodes)
		'style': {
			'opacity': 1,
			'label': label,
			'background-color': negative_color,
		}
	},
	{
		'selector': 'edge',						# for the edges
		'style': {
			'line-color': '#C5D3E2',			# very ligth blue
			'curve-style': 'haystack'
		}
	}
]

# Necessary to get more cytoscape layouts
cyto.load_extra_layouts()

"""
#####################
##  Middle window  ##
#####################
"""

# This will appear above the cytoscape node graph
node_graph_layout = dbc.Row(
	[
		dbc.Col(
			[
				# Drop down selector for layout of cytoscape, at the top of the page
				html.Main("Choose graph layout style:"),
				dcc.Dropdown(
					id='dropdown-layout',
					options=[
						{'label': 'random',
						'value': 'random'},
						{'label': 'grid',
						'value': 'grid'},
						{'label': 'circle',
						'value': 'circle'},
						{'label': 'concentric',
						'value': 'concentric'},
						{'label': 'breadthfirst',
						'value': 'breadthfirst'},
						{'label': 'cose',
						'value': 'cose'},
						{'label': 'euler',
						'value': 'euler'},
						{'label': 'cose-bilkent',
						'value': 'cose-bilkent'},
						{'label': 'cola',
						'value': 'cola'},
						{'label': 'spread',
						'value': 'spread'},
						{'label': 'dagre',
						'value': 'dagre'},
						{'label': 'klay',
						'value': 'klay'}
					], value='grid'
				),

				# Legend of the logFC colors in the cytoscape node graph
				dbc.Row([
					dbc.Col(html.Main("Log fold change: ")),
					dbc.Col([
						html.Div([
							html.Main("POSITIVE: "),
							html.Span(style={"height": "15px", "width": "15px", "border-radius": "50%", "display": "inline-block", "background-color": positive_color, "margin": "5px 10px"}),
						], className='row')
					]),
					dbc.Col([
						html.Div([
							html.Main("NEGATIVE: "),
							html.Span(style={"height": "15px","width": "15px","border-radius": "50%","display": "inline-block","background-color": negative_color, "margin": "5px 10px"}),
						], className='row')
					]),
					dbc.Col([
						html.Div([
							html.Main("SIMILAR: "),
							html.Span(style={"height": "15px","width": "15px","border-radius": "50%","display": "inline-block","background-color": neutral_color, "margin": "5px 10px"}),
						], className='row')
					]),
				])
			]
		)
	]
)

# Cystoscape graph itselfs
node_graph = dbc.Row(
	[
		dbc.Col(
			html.Div(children=[
				cyto.Cytoscape(
					id='cytoscape-protein',
					elements = cy_edges + cy_nodes,
					style = {
						'width': '100%',				# Take up 100% of the width of the space it has been assigned
						'height': '87vh'				# 87vh = 87% of the total screen height
					},
					stylesheet = default_stylesheet,  	# stylesheet for nodes and edges
					responsive = True					# Changes size cystocape graph is browser window changes size
				)
			])
		)
	]
)

# columns to display in interaction table
interaction_table_columns = ['node1', 'node1_uniprot','node1_kegg', 'node2', 'node2_uniprot', 'node2_kegg', 'interaction', 'combined_score']

interaction_table =  dbc.FormGroup([
	# Putting button on top, as sometimes you can't press this if it's at the bottom, and the table is very long
	dbc.Button(
		id='all_button',
		n_clicks=0,
		children='Show all proteins',
		color='info',
		block=True
	),
	html.Br(),
			
	dash_table.DataTable(
		id='interaction_table',
		columns=[{'name': i, 'id': i, 'deletable': False} for i in interaction_table_columns],
		data = nodeDf.to_dict('records'),
		page_size = 25,					# 25 rows
		filter_action='native',
		filter_query='',
		sort_action='native',
		sort_mode='multi',
		sort_by=[],
		style_cell={'textAlign': 'left', 'maxWidth': '350px', 'whiteSpace': 'normal'},
		style_cell_conditional=(
		[
			{
				'if': {'column_id': c},
				'backgroundColor': '#ffffe5'
			} for c in ['node1', 'node1_uniprot','node1_kegg']
		]+
		[
			{
				'if': {'column_id': c},
				'backgroundColor': '#ffffcc'
			} for c in ['node2', 'node2_uniprot','node2_kegg']
			
		]+
		[
			{
				'if': {'column_id': 'combined_score'},
				'textAlign': 'right',
			}
		]),
	),
])

# columns to display in acetylation table
acetylation_table_columns = ['geneName', 'uniprotID','numAcSites', 'peptides', 'protLogFC', 'keggPathways']

acetylation_table =  dbc.FormGroup([
	dbc.Button(
		id='all_button2',
		n_clicks=0,
		children='Show all proteins',
		color='info',
		block=True
	),
	html.Br(),
			
	dash_table.DataTable(
		id='acetylation_table',
		columns=[{'name': i, 'id': i, 'deletable': False} for i in acetylation_table_columns],
		
		#tooltip_header={i: i for i in acetylation_table_columns}
		# TODO: Column tooltips for extra information??

		data = acetylation.to_dict('records'),
		page_size = 25,
		filter_action='native',
		filter_query='',
		sort_action='native',
		sort_mode='multi',
		sort_by=[],
		style_cell={'textAlign': 'left', 'maxWidth': '350px', 'whiteSpace': 'normal'},
		style_cell_conditional=(
		[
			{
				'if': {'column_id': 'numAcSites'},
				'textAlign': 'center',
			}
		])
	),
	
])

# Operators for filtering data from DataTables
operators = [['ge ', '>='],
			 ['le ', '<='],
			 ['lt ', '<'],
			 ['gt ', '>'],
			 ['ne ', '!='],
			 ['eq ', '='],
			 ['contains ']]

# Used in callback for filtering DataTables
# returns column name where filter was applied, the operator used, and the filtervalue itself
def split_filter_part(filter_part):
	for operator_type in operators:
		for operator in operator_type:
			if operator in filter_part:
				name_part, value_part = filter_part.split(operator, 1)
				name = name_part[name_part.find('{') + 1: name_part.rfind('}')]
				value_part = value_part.strip()
				v0 = value_part[0]
				if (v0 == value_part[-1] and v0 in ("'", '"', '`')):
					value = value_part[1: -1].replace('\\' + v0, v0)
				else:
					try:
						value = float(value_part)
					except ValueError:
						value = value_part
				return name, operator_type[0].strip(), value
	return [None] * 3

# Middle window contains cytoscape graph and tables
middle_window = html.Div(
	id = 'middle-window',
	children=[
		# load all components that you want to display at some point!! Else you'll get errors
		node_graph_layout,
		node_graph,
		interaction_table,
		acetylation_table
	],
	style=
	{
	'margin-left': '16%',   # 16% from left side, as the sidebar takes up 15% + 1% for padding (looks nicer)
	'margin-top': '1%',
	'margin-right': '21%',
	'padding': '20px 10p'
	}
)



"""
#######################
##  Left side panel  ##
#######################
"""

# contains all components from left_side_panel that have a some function
# check callback function of these components for actual functionality
controls = dbc.FormGroup(
	[
		# Navigation buttons to switch between tables and cytoscape node graph
		dbc.Nav([
			dbc.NavLink("Interaction data table", href="/interaction_table", id="interaction_table-link"),
			dbc.NavLink("Acetylation data table", href="/acetylation_table", id="acetylation_table-link"),
			dbc.NavLink("Cytoscape (graphical view)", href="/cytoscape", id="cytoscape-link"),	
		],vertical=True),
		html.Hr(),
		# Node labels selector
		html.Div('Node labels:'),
		dcc.Dropdown(
			id='change_label',
			options=[
				{'label': 'Gene name',
				'value': 'pref_name'},
				{'label': 'StringDB',
				'value': 'StringDB'},
				{'label': 'Uniprot',
				'value': 'Uniprot'},
				{'label': 'KEGG ID',
				'value': 'KEGG ID'},
			], value='pref_name'
		),
		# Some space between components
		html.Br(),

		# search function input bar and button
		html.Div('Search a protein:'),
		dcc.Input(
			id='searchvalue',
			type='search',
			placeholder='Search protein',
			style={'maxWidth':'50%'}
		),
		html.Button(id='searchbutton', n_clicks=0, children='Search'),

		html.Br(),
		html.Hr(),

		# Callback for this will display the interactions attribute for edges as an edge label
		dcc.Checklist(
			id= 'edgelabel-options',
			options=[
				{'label': '  Show interaction label on edges', 'value': 'show_interaction'}
			]
		),

		html.Br(),
		
		# Button that hides all nodes that do not have any annotated data // shows them again when pressed a second time
		html.Div('Switch between full or subset of data:'),
		dbc.Button(
			id='unique_button',
			n_clicks=0,
			children='Show only annotated proteins',
			color='secondary',
			block=True
		),
		html.Br(),

		dbc.Button("Return to default look", id="to_default_stylesheet", block=True, color='primary'),

		# Exports cytoscape node graph as an image
        html.Div([
			html.Main('Export graph as image:'),
			html.Div([
				dbc.Button("JPG", id="btn-get-jpg",block=False,color='primary', className="mr-2"),
            	dbc.Button("PNG", id="btn-get-png",block=False,color='info', className="mr-2"),
            	dbc.Button("SVG", id="btn-get-svg",block=False,color='primary', className="mr-2")
			])
        ], 
		className='col',
		style= {
			'position': 'fixed',
			'bottom': 10,
			'padding': '20pw 10px'
		}),

    ]
)

# Left side panel that contains all the controls
left_side_panel = html.Div(
	[
		html.H3(
			'PACES', 
			style= 
			{
				'textAlign': 'center',
				'color': '#191970'
			}
		),
		html.Hr(),	# Draws horizontal line
		controls,
	],
	style=
		{
		'position': 'fixed',
		# positioning
		'top': 0,
		'left': 0,
		'bottom': 0,
		# width of the sidebar == 15% of the browser window's width
		'width': '15%',
		# some padding to make it look nicer
		'padding': '20px 10px',
		'background-color': '#f8f9fa'
		}
)

# Defined some default text styles
CARD_TEXT_STYLE = {
	'margin-left': 10,
	'margin-right': 5
}
CARD_TEXT_STYLE_UNDERLINE = {
	'margin-left': 10,
	'margin-right': 5,
	'text-decoration': 'underline'
}


"""
########################
##  Right side panel  ##
########################
"""

# Right side panel contains all information for the selected node from the cytoscape graph
right_side_panel = html.Div(
    [
        dbc.FormGroup(
            [
				
				dbc.Card(
					children=[
						dbc.CardHeader([html.H4('Selected node details')]),
						dbc.CardBody([

							# Protein ID
							html.Div([
								html.Main('Gene name: ', style=CARD_TEXT_STYLE),
								html.Main('Nothing selected', id = 'selectedNode-id', style={'margin-left': 1})
							], className='row'),

							# Uniprot ID
							html.Div([
								html.Main('Uniprot ID: ', style=CARD_TEXT_STYLE),
								html.A('Nothing selected', href='', id='selectedNode-uniprot', target="_blank", style={'margin-left': 7})	# Link to database
							], className='row'),

							# STRING ID
							html.Div([
								html.Main('STRING ID: ', style=CARD_TEXT_STYLE),
								html.A('Nothing selected', href='', id = 'selectedNode-string', target="_blank", style={'margin-left': 7})	# Link to database
							], className='row'),

							# KEGG ID
							html.Div([
								html.Main('KEGG ID: ', style=CARD_TEXT_STYLE),
								html.A('Nothing selected', href='', id = 'selectedNode-kegg', target="_blank", style={'margin-left': 21})	# Link to database
							], className='row'),

							html.Div([
								html.Main('Acetylated sites: ', style=CARD_TEXT_STYLE),
								html.Main('Nothing selected', id = 'selectedNode-acetylation_sites', style={'margin-left': 1})
							], className='row'),

							html.Div([
								html.Main('logFC: ', style=CARD_TEXT_STYLE),
								html.Main('Nothing selected', id = 'selectedNode-logFC', style={'margin-left': 5})
							], className='row'),

							html.Div([
								html.Hr()
							], className='row'),

							# Annotation data
							html.Div([
								html.Main('Annotation data:', style=CARD_TEXT_STYLE_UNDERLINE)
							], className='row'),

							html.Div([
								html.Main('Nothing selected', id='selectedNode-annotation', style=CARD_TEXT_STYLE)
							], className='row'),

							html.Div([
								html.Hr()
							], className='row'),

							# Annotation data
							html.Div([
								html.Main('KEGG pathways:', style=CARD_TEXT_STYLE_UNDERLINE)
							], className='row'),


                            html.Div([
                                dcc.Markdown('Nothing selected', id='selectedNode-kegganno', style=CARD_TEXT_STYLE, dangerously_allow_html=True)
                            ], className='row'),
							
							html.Main("Currently displaying {} nodes ".format(len(set(nodes))), id='total_nodes')
                        ],
                        style={'max-height': '85vh', 'overflow-y': 'auto'},	# overflow = auto --> makes card scrollable if text is to large to fit screen
                        )
                    ]
                ),
            ]
        )
    ],
    style =
		{
		'position': 'fixed',
		'top': 0,
		'right': 0,
		'bottom': 0,
		'width': '20%',
		'padding': '20px 10px',
		'background-color': '#f8f9fa'
		}
)


"""
############################
## Define Dash app itself ##
############################
"""

# Define app itself
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])

# Set app layout
app.layout = html.Div([dcc.Location(id="url"), left_side_panel, right_side_panel, middle_window])


"""
####################
## Dash Callbacks ##
####################
"""

# Changes layout of the cytoscape node graph
@app.callback(
	Output('cytoscape-protein', 'layout'),
	[Input('dropdown-layout', 'value')
])
def update_cytoscape_layout(layout):
	return {'name': layout}


# Displays total number of nodes in cytoscape node graph
@app.callback(Output('total_nodes','children'),
			Input('cytoscape-protein', 'elements'))
def num_nodes(elz):
	nodes=[]
	for item in elz:
		if 'source' in str(item):
			nodes.append(item['data']['source'])
			nodes.append(item['data']['target'])
	return "Currently displaying {} nodes ".format(len(set(nodes)))


# This function gets all information of the selected node, converts it to strings, makes links to databases, and displays it in the right side panel
@app.callback( Output('selectedNode-id', 'children'),
				Output('selectedNode-uniprot', 'children'),
				Output('selectedNode-uniprot', 'href'),
				Output('selectedNode-string', 'children'),
				Output('selectedNode-string', 'href'),
				Output('selectedNode-kegg', 'children'),
				Output('selectedNode-kegg', 'href'),
				Output('selectedNode-acetylation_sites', 'children'),
				Output('selectedNode-logFC', 'children'),
				Output('selectedNode-annotation', 'children'),
				Output('selectedNode-kegganno', 'children'),
				[Input('cytoscape-protein', 'tapNodeData')])
def displaySelectedNodeData(data):   
# data is a dictionary, can be returned as json with 'return json.dumps(data, indent=2)'

	if data:	# This is necessary, if not here, then data will not be dictionary but a NoneType
		prot_id = str(data.get('id'))
		uniprot_id = str(data.get('UniprotID'))
		uniprot_link = "https://www.uniprot.org/uniprot/{}".format(uniprot_id)
		string_id = str(data.get('stringid'))
		string_link = "https://string-db.org/network/{}".format(string_id)
		kegg_id = str(data.get('KEGG_ID'))
		kegg_link = "https://www.genome.jp/dbget-bin/www_bget?pae:{}".format(kegg_id)

		acetylation_sites = str(acet_sites.get(uniprot_id))
		logfc = str(data.get('logFC'))

		if prot_id in unique_keys:
			annotation = str(protein_annotation.get(string_id))
		else:
			annotation = "Annotation not available"
		
		# (Try, expect) has to be used here else a TypeError will occur when there is no KEGG pathways for the node
		try:
			kegg_annotation = '<br>'.join(kegg_dict.get(kegg_id))
		except TypeError:
			kegg_annotation = "No path"

		return prot_id, uniprot_id, uniprot_link, string_id, string_link, kegg_id, kegg_link, acetylation_sites, logfc, annotation, kegg_annotation
	else:
		return 'Nothing selected', 'Nothing selected', '', 'Nothing selected', '', 'Nothing selected', '', 'Nothing selected', 'Nothing selected', 'Nothing selected', 'Nothing selected'				


# Makes new cy_edges cy_nodes with only nodes that have annotation data to pass to cytoscape graph
@app.callback(
	Output('cytoscape-protein', 'elements'),
	Output('unique_button', 'children'),
	[Input('unique_button', 'n_clicks')
])	
def only_show_annotated_cytoscape(clix):
	ctx = dash.callback_context
	nodes=set()
	cy_edges = []
	cy_nodes = []

	if ctx.inputs['unique_button.n_clicks']%2 == 1:
		#If its the unique nodes		
		for index, row in nodeDf.iterrows():
			source, target, score, interaction = row['node1'], row['node2'], row['combined_score'], row['interaction']
			source_uniprot = row['node1_uniprot']
			target_uniprot = row['node2_uniprot']

			source_logFC = get_logFC_as_string(logfc_anno.get(source_uniprot))
			target_logFC = get_logFC_as_string(logfc_anno.get(target_uniprot))

			if source not in nodes and source in unique_keys:
				nodes.add(source)
				cy_nodes.append({'data': {'id': source, 'label': source, 'UniprotID': source_uniprot, 'logFC': source_logFC}})	  # data: { [here comes that data contained for the node as a dictionary] }
			if target not in nodes and target in unique_keys:
				nodes.add(target)
				cy_nodes.append({'data': {'id': target, 'label': target, 'UniprotID': target_uniprot, 'logFC': target_logFC}})

			if row['combined_score'] >= NODE_CUTOFF_SCORE and source in unique_keys and target in unique_keys:
				cy_edges.append({'data': {'source': source, 'target': target, 'score': score, 'interaction': interaction, 'source_logFC': source_logFC, 'target_logFC': target_logFC}})	 # same for the edge	
			
	else:
			#Return all nodes back to original graph
		for index, row in nodeDf.iterrows():
			source, target, score, interaction = row['node1'], row['node2'], row['combined_score'], row['interaction']
			source_uniprot = row['node1_uniprot']
			target_uniprot = row['node2_uniprot']
			
			source_logFC = get_logFC_as_string(logfc_anno.get(source_uniprot))
			target_logFC = get_logFC_as_string(logfc_anno.get(target_uniprot))

			if source not in nodes:
				nodes.add(source)
				cy_nodes.append({'data': {'id': source, 'label': source, 'UniprotID': source_uniprot, 'logFC': source_logFC}})	  # data: { [here comes that data contained for the node as a dictionary] }
			if target not in nodes:
				nodes.add(target)
				cy_nodes.append({'data': {'id': target, 'label': target, 'UniprotID': target_uniprot, 'logFC': target_logFC}})

			cy_edges.append({'data': {'id': source+target, 'source': source, 'target': target, 'score': score, 'interaction': interaction, 'source_logFC': source_logFC, 'target_logFC': target_logFC}})	 # same for the edge	

	button_text = 'Show only annotated proteins' if ctx.inputs['unique_button.n_clicks']%2 == 0 else 'Show all proteins'
	
	elements = cy_edges + cy_nodes
	return elements, button_text
				

# Changes layout of middle-window depending on url
@app.callback(Output('middle-window', 'children'),
				[Input('url', 'pathname')])
def updateMiddleWindow(pathname):
	if pathname == '/':
		page=html.Div([html.P(
		html.H3("Welcome to PACES webpage by Adi, Ben, Stefaan & Hannelore!"),style={'color':'#007fff','text-align':'center'}),
		html.Br(),
		html.P("Note: Any changes made to the tabular view (Interaction data table/ Acetylation data table) automatically reflect on the graphical view (Cytoscape)",style={'text-align':'center'})
		])
		return page
	if pathname == '/cytoscape':
		return (node_graph_layout, node_graph)
	if pathname == '/interaction_table':
		return interaction_table
	if pathname == '/acetylation_table':
		return acetylation_table
	
	
# callback for interaction table
@app.callback(
	Output('interaction_table', 'data'),
	[Input('interaction_table', 'sort_by'),
	 Input('interaction_table', 'filter_query'),
	 Input("all_button", "n_clicks")
	 ])
def update_interaction_table(sort_by, filter,n_clicks):
	global nodeDf
	global acetylation

	changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
	if 'all_button' in changed_id:
		nodeDf = nodeDf_orig
		acetylation = acetylation_orig
		
	filtering_expressions = filter.split(' && ')

	dff = nodeDf
	# Do this before search, else some results might get excluded
	dff = dff.round({'combined_score': 3})
	for filter_part in filtering_expressions:
		col_name, operator, filter_value = split_filter_part(filter_part)
		if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
			# these operators match pandas series operator method names
			dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
		elif operator == 'contains':	
			if 'kegg' in col_name:
				dff_na = dff.dropna()
				dff = dff_na.loc[dff_na[col_name].str.contains(str(filter_value))]
			else:
				dff = dff.loc[dff[col_name].str.contains(str(filter_value))]
			
	if len(sort_by):
		dff = dff.sort_values(
			[col['column_id'] for col in sort_by],
			ascending=[
				col['direction'] == 'asc'
				for col in sort_by
			],
			inplace=False
		)

	# make the filtered dataframe the nodeDf which is then also used in the cytoscape graph
	nodeDf=dff

	# updates acetylation dataframe with filtered entries from nodeDf (so filter from interatction table also applies to acetylation table)
	proteins = sum([ nodeDf['node1_uniprot'].tolist(), nodeDf['node2_uniprot'].tolist()], [])
	acetylation = acetylation[acetylation['uniprotID'].isin(proteins)]

	# Rounds float of combined score in table view
	dff = dff.round({'combined_score': 3})

	return dff.to_dict('records')

# callback for acetylation table
@app.callback(
	Output('acetylation_table', 'data'),
	[Input('acetylation_table', 'sort_by'),
	 Input('acetylation_table', 'filter_query'),
	 Input("all_button2", "n_clicks")
	 ])
def update_acetylation_table(sort_by, filter,n_clicks):
	global acetylation
	global nodeDf
		
	filtering_expressions = filter.split(' && ')

	changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
	if 'all_button2' in changed_id:
		nodeDf = nodeDf_orig
		acetylation = acetylation_orig

	dff = acetylation
	for filter_part in filtering_expressions:
		col_name, operator, filter_value = split_filter_part(filter_part)
		if operator in ('eq', 'ne', 'lt', 'le', 'gt', 'ge'):
			# these operators match pandas series operator method names
			dff = dff.loc[getattr(dff[col_name], operator)(filter_value)]
		elif operator == 'contains':
			dff = dff.loc[dff[col_name].str.contains(str(filter_value))]
			
	if len(sort_by):
		dff = dff.sort_values(
			[col['column_id'] for col in sort_by],
			ascending=[
				col['direction'] == 'asc'
				for col in sort_by
			],
			inplace=False
		)
	# make the filtered dataframe the acetlylation dataframe which is then also used to filter the nodeDf
	acetylation=dff
	
	# updates nodeDf with filtered entries from acetylation
	proteins = acetylation['uniprotID'].tolist()
	nodeDf = nodeDf[nodeDf['node1_uniprot'].isin(proteins) | nodeDf['node2_uniprot'].isin(proteins)]

	# Rounds float of combined score in table view
	dff = dff.round({'combined_score': 3})

	return dff.to_dict('records')

# Export current cytoscape graph as an image when button is clicked
@app.callback(
	Output("cytoscape-protein", "generateImage"),
	[
		Input("btn-get-jpg", "n_clicks"),
		Input("btn-get-png", "n_clicks"),
		Input("btn-get-svg", "n_clicks"),
	])
def getImage(jpg, png, svg):
	ftype=''
	action = 'download'
	if dash.callback_context.triggered:
		input_id = dash.callback_context.triggered[0]['prop_id'].split('.')[0]
		if input_id != "tabs":
			action = "download"
			ftype = input_id.split("-")[-1]
	return {'type': ftype, 'action':action}


""" Change layout of cytoscape node graph

inputs used to change layout:
 - selecting a node:											Colors selected node, and the nodes and edges connected to it
 - clicking 'Return to default look' button:					Changes cytoscape stylesheet to the default look
 - selecting a different label in the node labels dropdown:		Changes label of nodes
 - clicking search button:										Uses state of the 'searchvalue' field to look for a node and color it purple
 - clicking edgelabel checklist:								Displays interactions atributes of an edge as an edge label.
"""
@app.callback(Output('cytoscape-protein', 'stylesheet'),
			  [Input('cytoscape-protein', 'tapNode'),
			  Input("to_default_stylesheet", "n_clicks"),
			  Input('change_label', 'value'),
			  Input('searchbutton', 'n_clicks'),
			  Input('edgelabel-options', 'value'),
			  State('searchvalue', 'value'),])
def generate_stylesheet(node, button, new_label, searchbutton, edgelabelvalue, searchvalue):

	global label						# globally changes label value (= also outside of this function)
	if new_label == 'pref_name':
		label = 'data(label)'
	if new_label == 'StringDB':
		label = 'data(stringid)'
	if new_label == 'Uniprot':
		label = 'data(UniprotID)'
	if new_label == 'KEGG ID':
		label = 'data(KEGG_ID)'
	
	# Default stylesheet has to be defined again so that label value is updated
	new_default_stylesheet = default_stylesheet
	
	# returns to default stylesheet
	changed_id = [p['prop_id'] for p in dash.callback_context.triggered][0]
	if 'to_default_stylesheet' in changed_id:
		return new_default_stylesheet

	# if no node is selected, and search button is not selected, return to default stylesheet
	if not node and not 'searchbutton' in changed_id:
		return new_default_stylesheet

	# Defines value of edge label
	edge_label = ''
	try:
		if 'show_interaction' in edgelabelvalue:
			edge_label = "data(interaction)"	
	except:
		pass
	
	# beginning of new stylesheet. Everything is very opaque
	stylesheet = [{
		"selector": 'node',
		'style': {
			'opacity': 0.3
		}
	}, {
		'selector': 'edge',
		'style': {
			'opacity': 0.2,
			"curve-style": "bezier",
		}
	}]

	# If searchbutton was pressed, change color of node that matches search value
	if 'searchbutton' in changed_id:
		search_id = ''
		if str(searchvalue) in nodes:
			search_id = 'id'
		if str(searchvalue)[0].isalpha() and str(searchvalue[1]).isdigit():
			search_id = 'UniprotID'
		if 'DR' in str(searchvalue).upper():
			search_id = 'stringid'
		if 'PA' in str(searchvalue).upper():
			search_id = 'KEGG_ID'

		if not search_id:
			return stylesheet
		else:
			stylesheet.append({
					"selector": 'node[{} = "{}"]'.format(search_id, str(searchvalue)),
					"style": {
						'background-color': '#B10DC9',
						"border-color": "purple",
						"border-width": 2,
						"border-opacity": 1,
						"opacity": 1,

						"label": label,
						"text-opacity": 1,
						'z-index': 9999			# the bigger the Z-index, the more priority it has to be in front of other nodes/labels
					}
				})
			return stylesheet

	# if searchbutton was not pressed
	elif 'searchbutton' not in changed_id:
		stylesheet.append({
			"selector": 'node[id = "{}"]'.format(node['data']['id']),
			"style": {
				'background-color': '{}'.format(colordict.get(node['data']['logFC'])),
				"border-color": "purple",
				"border-width": 2,
				"border-opacity": 1,
				"opacity": 1,

				"label": label,
				"text-opacity": 1,
				'z-index': 9999
			}
		})
		

	for edge in node['edgesData']:
		if edge['source'] == node['data']['id']:
			stylesheet.append({
				"selector": 'node[id = "{}"][logFC = "{}"]'.format(edge['target'], edge['target_logFC']),
				"style": {
					'background-color': '{}'.format(colordict.get(edge['target_logFC'])),
					"opacity": 1,
					"label": label,
					"text-opacity": 1,
					'z-index': 9999
				}
			})
			stylesheet.append({
				"selector": 'edge[id= "{}"]'.format(edge['id']),
				"style": {
					"label": "{}".format(edge_label),
					"text-rotation": "autorotate",
					"line-color": selected_edge_color,
					'opacity': 1,
					'z-index': 5000
				}
			})

		if edge['target'] == node['data']['id']:
			stylesheet.append({
				"selector": 'node[id = "{}"][logFC = "{}"]'.format(edge['source'], edge['source_logFC']),
				"style": {
					'background-color': '{}'.format(colordict.get(edge['source_logFC'])),
					"opacity": 1,

					"label": label,
					"text-opacity": 1,
					'z-index': 9999
				}
			})
			stylesheet.append({
				"selector": 'edge[id= "{}"]'.format(edge['id']),
				"style": {
					"label": "{}".format(edge_label),
					"text-rotation": "autorotate",
					"line-color": selected_edge_color,
					'opacity': 1,
					'z-index': 5000
				}
			})
	return stylesheet
		

# run app
if __name__ == '__main__':
	app.run_server(debug=True)