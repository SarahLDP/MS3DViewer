# Standard library imports
import socket
import time

import subprocess
import sys

import base64
import os



def install_requirements(filename='requirements.txt'):
    """
    Install dependencies listed in a requirements text file.

    Parameters:
    - filename (str): The name of the requirements text file. Default is 'requirements.txt'.

    Raises:
    - CalledProcessError: If the installation process fails.

    Example:
    install_requirements('requirements.txt')
    """
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    requirements_path = os.path.join(script_dir, filename)

    subprocess.check_call([sys.executable, "-m", "pip", "install", "-r", requirements_path])

install_requirements()

# Third-party library imports

from dash import Dash, dcc, html, ctx
from dash.dependencies import Input, Output, State
from tqdm import tqdm
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.pyplot as plt

# Custom module imports
import parse_file
from classes import Protein, Peptide
from viewer import create_viewer
import peptide_atlas
import get_colors
import info

# ######################################################
    ### Create class objects ###

protein_list = []
peptide_list = []
file_path = ''

def create_class_objs(workbook,protein_df,protein_index):
    """
    Creates instances of Protein and Peptide classes based on data from a workbook.

    Parameters:
    - workbook (Workbook): The workbook containing protein and peptide data.
    - protein_df (DataFrame): DataFrame containing protein information.
    - protein_index (list of int): List containing the starting index of each protein's data in the DataFrame,
      used to locate and extract relevant peptide data.

    Returns:
    - tuple: A tuple containing:
        - list: List of Protein objects.
        - list: List of Peptide objects.

    Note:
    - The Protein and Peptide classes must be defined. 
    """
    global protein_list
    global peptide_list

    i=0

    for index, row in protein_df.iterrows():
        i+=1

        protein = Protein(
            accession=row['Accession'],
            total_psms=row['# PSMs']     
        )
        protein_list.append(protein)

        print(f'Protein {protein.accession} {i} of {len(protein_df)}')


        protein_peptides = []
        peptide_df = parse_file.extract_peptide_df(workbook,protein.accession,protein_df,protein_index)
        # Wrap the pandas iterrows with tqdm for the progress bar
        for pep_index, pep_row in tqdm(peptide_df.iterrows(), total=peptide_df.shape[0], desc=f'Processing peptides for protein {protein.accession}'):
            moddict = {}
            for r in pep_row['Modifications']:
                moddict[r[0]] = r[1:]  # Creating a dictionary of modifications for each peptide

            peptide = Peptide(protein=protein,  # Creating peptide objects of the class Protein
                sequence=pep_row['Annotated Sequence'],
                start_position=int(pep_row['Positions in Proteins'][0]), 
                end_position=int(pep_row['Positions in Proteins'][1]),
                modifications=moddict,
                file_id = pep_row['File ID']
                )
            peptide_list.append(peptide)
            protein_peptides.append(peptide)
        protein.set_peptides(protein_peptides)
    return protein_list, peptide_list



######################################

### Creating the Dash app ###

app = Dash(__name__)


def find_available_port(start_port=5050, max_tries=10):
    """
    Finds an available port within a specified range.

    Parameters:
    - start_port (int): The starting port number to begin the search. Default is 5050.
    - max_tries (int): The maximum number of port attempts to make. Default is 10.

    Returns:
    - int: The available port number if found.

    Raises:
    - OSError: If no available port is found within the specified range.

    Example:
    >>> port = find_available_port()
    >>> print(port)

    Note:
    - This function attempts to bind to each port in the specified range until an available port is found.
    - It raises an OSError if no available port is found within the specified number of tries.
    """
    for port in range(start_port, start_port + max_tries):
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.bind(("localhost", port))
            return port
        except OSError:
            pass
    raise OSError("Could not find an available port")


def view_3d(v):
    """
    Generates an HTML page containing a 3D viewer object.

    Parameters:
    - v: The 3D viewer object.

    Returns:
    - str: An HTML page containing the 3D viewer object.

    Example:
    >>> html_page = view_3d(viewer_object)
    >>> print(html_page)

    Note:
    - This function generates an HTML page embedding a 3D viewer object using the provided viewer object.
    - It assigns a unique identifier to the viewer object for JavaScript interaction.
    - JavaScript functions are embedded in the HTML page to interact with the viewer object.
    """
    v.uniqueid = str(time.time_ns())
    js_functions = """
    <script>
        function getViewerState() {
            return viewer_""" + v.uniqueid + """.getView();
        }
    </script>
    """
    return '''<html style="overflow: hidden; height: 100%; width: 100%"> \
                <head>
            
            <script src="http://code.jquery.com/jquery-2.1.3.min.js"></script>
          
            <script>
             "use strict";
             let myview
             $(window).load(() => {
                  document.getElementById("undefined").addEventListener("click", () => console.log( viewer_''' + v.uniqueid + '''.getView()))
                 })
            </script>
        </head>
     <body style="height: 100%; width: 100%">''' + v._make_html() + js_functions + '</body></html>'




app.layout = html.Div(children=[
    html.Div(id='upload-section',children=[
        dcc.Upload(
            id='upload-data',
            children=html.Button('Upload File'),
            multiple=False
        )
    ]),
    html.Div(id='output-data-upload'),
    html.Div(id='main-layout'),
    html.Hr()
])


@app.callback(
    Output('output-data-upload', 'children'),
    Input('upload-data', 'contents'),
    State('upload-data', 'filename')
)
def save_upload(contents,filename):
    """
    Process the uploaded file contents and save it to a local directory.
    
    Args:
        contents (str): Base64 encoded contents of the file.
        filename (str): Name of the file to be saved.

    Returns:
        A layout component showing the main application layout if a file is uploaded,
        otherwise a message indicating no file has been uploaded.
    """
    global file_path

    if contents is not None:
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string)
        file_path = os.path.join('uploaded_files', filename)

        if not os.path.exists('uploaded_files'):
            os.makedirs('uploaded_files')

        with open(file_path, 'wb') as f:
            f.write(decoded)
        
        return html.Div([f"""File uploaded: {file_path}."""])
    return html.Div(['No file uploaded.'])


@app.callback(
        Output('main-layout','children'),
        Output('upload-section','children'),
        Input('output-data-upload', 'children')
)
def get_main_layout(content):
    if not content==html.Div(['No file uploaded.']) and not file_path=='':
        return main_app_layout(file_path),html.Div('File Uploaded')



def main_app_layout(file_path):
    """
    Generate the main application layout based on the contents of the uploaded file.
    
    Args:
        file_path (str): The path to the uploaded file processed for visualization.

    Returns:
        An HTML component containing the visual layout elements for protein and peptide data visualization.
    """
    ## workbook ##
    global protein_list
    global peptide_list
    global samples_in_file
    workbook = parse_file.open_workbook(file_path)
    protein_df, protein_index, accession_numbers = info.get_protein(workbook)


    samples_in_file = info.get_samples_in_file(workbook)

    protein_list,peptide_list= create_class_objs(workbook=workbook,protein_df=protein_df,protein_index=protein_index)


    mod_in_fst_protein = []
    for m in protein_list[0].get_protein_modification_types(nterm=False).values():
        for mod in m:
            mod_in_fst_protein.append(mod)

    return html.Div([
        html.Div(id='output-data-upload'),
        html.Hr(),
        html.H3('Choose sample'),
        dcc.Dropdown(['All Samples']+samples_in_file, 'All Samples',id = 'sample-dropdown'),
        
        html.H3('Choose accession to display'),
        dcc.Dropdown(accession_numbers, accession_numbers[0],id = 'accession-dropdown'),
        
        html.Div([
            html.Div(children=[
                html.Div(id='modifications-container',
                    children=[
                        html.H3('Choose modifications: '),
                        dcc.Checklist(
                            options=protein_list[0].get_protein_modification_types(nterm=False),id='modification-input'
                            ),
                        html.Button('Update Viewer', id='apply-modifications-btn',style = {'font-family': 'Georgia','font-size':'16px'})
                    ], style={'padding': 10}),
                html.Br(),
                html.Label('Show residue labels: '),
                dcc.RadioItems(options=[{'label': 'On', 'value': 'On'}, {'label': 'Off', 'value': 'Off'}],
                            value = 'On',
                            id='label-choice'),
                html.Br(),
                dcc.Checklist(options = ['Show N-Term Modifications'],id = 'show-nterm-check'),
                dcc.Checklist(options = ['Remove Leader MET'],value=['Remove Leader MET'],id = 'remove-met-check'),
                # html.Label(id='met-test-label'),
                html.Br(),html.Br(),
                html.Label('Choose backbone style'),
                dcc.Dropdown(['Cartoon','Stick','Sphere','Cross','Line'],'Cartoon',id = 'vis_style_selector'),
                html.Br(),
                html.Label('Choose backbone color'),
                dcc.Dropdown(['By peptide abundance','Rainbow','Grey'],'By peptide abundance',id = 'vis_color_selector'),
                html.Br(),
                html.Label('Choose backbone size'),
                dcc.Input(id = 'vis_size_selector',type='number',value=0.4,min=0.1,max=10,step=0.1,style = {'font-family': 'Times New Roman','font-size':'15px'}),
                html.Br(),html.Br(),
    
                html.Label('Choose residue style'),
                dcc.Dropdown(['Cartoon','Stick','Sphere','Cross','Line'],'Stick',id = 'res_style_selector'),
                html.Br(),
                html.Label('Choose residue size'),
                dcc.Input(id = 'res_size_selector',type='number',value=0.4,min=0.1,max=10,step=0.1,style = {'font-family': 'Times New Roman','font-size':'15px'}),
                html.Br()
                ],style={'width': '15%', 'height': '650px'}),
            
           
            html.Div(children=[
                html.Iframe(
                    # srcDoc = view_3d(create_viewer(protein_list[0],protein_list[0].peptides,nterm=False)),
                    id='molView', 
                    style={'width': '100%', 'height': '650px'}),
                html.Div(id='color_bars')
            ],style={'width': '64%'}),
            
            html.Div(id='view_text',children=[],style = {'width':'20%','height':'650px'})
            
            ],style={'display':'flex','flexDirection':'row'}),
        html.Br(),   
        html.Label(id='Click-info'),
        html.Br(), html.Br(),  
        html.Button('Reset zoom',id='Zoom-reset',style = {'font-family': 'Georgia','font-size':'16px'}),
        html.Div(id='Click-data',style={'display':'none'}),
        html.Div([
            dcc.Graph(figure= peptide_atlas.peptide_atlas(protein_list[0].peptides,nterm=False)[0],id='pepView', style={'width': '100%', 'height': '700px'})
        ])
    ])





@app.callback(
    Output('pepView','figure'),
    Input('accession-dropdown','value'),
    Input('remove-met-check','value'),
    Input('show-nterm-check','value'),
    Input('sample-dropdown','value')
)
def update_peptide_atlas(accession,remove_m1,nterm,sample):
    """
    Update the peptide atlas visualization based on the selected accession number, sample, and modification options.

    Args:
        accession (str): The accession number of the protein to display.
        remove_m1 (bool): Whether to remove the initial methionine (M1) from the sequence visualization.
        nterm (bool): Whether to show N-terminal modifications in the visualization.
        sample (str): The sample file ID used to filter the peptides. If 'All Samples' is selected, all peptides are shown.

    Returns:
        A plotly graph object representing the updated peptide atlas visualization.
    """
    
    for protein in protein_list:
        if accession == protein.accession:
            s_protein = protein
    
    if sample == 'All Samples':
        peptides = s_protein.peptides

    else:
        peptides = s_protein.get_peptides_by_file_id(sample)
            
    range_max = len(s_protein.master_sequence) + 1
    
    if remove_m1:
        pa, pg, grouped = peptide_atlas.peptide_atlas(peptides=peptides, gap=1, range_max=range_max,remove_m1=True,nterm=nterm)
        pa = peptide_atlas.add_modification_legend(pa, grouped, s_protein, remove_m1=True,nterm=nterm)
    else:
        pa, pg, grouped= peptide_atlas.peptide_atlas(peptides=peptides, gap=1, range_max=range_max,remove_m1=False,nterm=nterm)
        pa = peptide_atlas.add_modification_legend(pa, grouped, s_protein, remove_m1=False,nterm=nterm)
    
    return pa


@app.callback(
        Output('modifications-container','children'),
        Input('accession-dropdown','value'),
        Input('show-nterm-check','value')
)
def update_modifications_list(accession,nterm):
    """
    Dynamically updates the list of protein modifications available for the selected protein accession number,
    based on whether N-terminal modifications are to be shown.

    Args:
        accession (str): The accession number of the protein for which modifications are to be listed.
        nterm (bool): Flag indicating whether N-terminal modifications should be included.

    Returns:
        A Dash HTML component containing a checklist of modifications and an update button.
        If no accession is selected, returns an empty list.
    """
    if accession is None:
        return []  # If no accession is selected, return an empty list
    else:
        # Use the selected accession to get modifications options for that protein
        for protein in protein_list:
            if accession == protein.accession:
                moddict_values = [value for value in protein.get_protein_modification_types(nterm=nterm).values()]
                for v in moddict_values:
                    modifications_in_protein = v
        return html.Div(children=[
            html.H3('Choose modifications: '),
            dcc.Checklist(
                options=modifications_in_protein,
                value = [option for option in modifications_in_protein],
                id='modification-input'
            ),
            html.Br(),
            html.Button('Update Viewer', id='apply-modifications-btn',style = {'font-family': 'Georgia','font-size':'16px'})
        ])  



    # modifications_in_protein = []  # Ensure always initialized
    # if accession:
    #     for protein in protein_list:
    #         if protein.accession == accession:
    #             # Assuming get_protein_modification_types returns a dictionary of modification types
    #             mod_types = protein.get_protein_modification_types(nterm=nterm)
    #             modifications_in_protein = [value for mod_list in mod_types.values() for value in mod_list]
    #             break  # Exit after finding the first match

    # if not modifications_in_protein:
    #     return html.Div(['No modifications available for the selected protein.'])

    # return html.Div([
    #     html.H3('Choose modifications: '),
    #     dcc.Checklist(
    #         options=modifications_in_protein,
    #         value=[option for option in modifications_in_protein],
    #         id='modification-input'
    #     ),
    #     html.Br(),
    #     html.Button('Update Viewer', id='apply-modifications-btn', style={'font-family': 'Georgia', 'font-size': '16px'})
    # ]) 

@app.callback(
    Output('molView', 'srcDoc'),
    Output('view_text','children'),
    Output('color_bars','children'),
    Input('apply-modifications-btn', 'n_clicks'),
    State('modification-input', 'value'),
    Input('accession-dropdown','value'),
    Input('label-choice','value'),
    Input('remove-met-check','value'),
    Input('vis_style_selector','value'),
    Input('vis_color_selector','value'),
    Input('vis_size_selector','value'),
    Input('res_style_selector','value'),
    Input('res_size_selector','value'),
    Input('Click-data','children'),
    Input('sample-dropdown','value'),
    Input('show-nterm-check','value')
)
def update_view(n_clicks,selected_modifications,accession,label_choice,remove_m1,visstyle_value,viscolor_value,vissize_value,resstyle_value,ressize_value,peptide_click,sample,nterm,mod_freq_color = 'YlOrRd'):  
    """
    Updates the molecular viewer and information tabs based on the user's selections including protein modifications,
    visual styles, and other visualization settings.

    Args:
        n_clicks (int): Number of times the 'apply modifications' button was clicked.
        selected_modifications (list): List of selected modifications to display.
        accession (str): Accession number of the protein.
        label_choice (str): User's choice for labeling the residues.
        remove_m1 (bool): Whether to remove the first methionine from the visualization.
        visstyle_value (str): Style of visualization for the molecule (e.g., 'cartoon', 'stick').
        viscolor_value (str): Color scheme for the visualization.
        vissize_value (float): Size of the backbone in the visualization.
        resstyle_value (str): Style of visualization for the residues.
        ressize_value (float): Size of the residues in the visualization.
        peptide_click (list): Range of residues clicked for zooming.
        sample (str): Sample selection to filter the data.
        nterm (bool): Whether to include N-terminal modifications in the visualization.
        mod_freq_color (str, optional): Color scheme for frequency of modifications (default 'YlOrRd').

    Returns:
        tuple: Contains the updated source document for the molecular viewer and a div with visualization information.
    """
    visstyle_value = visstyle_value.lower()
    resstyle_value = resstyle_value.lower()
    
    # For each protein in the protein list
    for protein in protein_list:
        protein:Protein
        ## if the accession number of the protein is the chosen accession number, 
            ### check if the selected modifications are empty or not
        if protein.accession == accession:

            if sample == 'All Samples':
                pfound = "## Showing all samples"
                peptidelist = protein.peptides

                
                max_psm = 0
                psm_mod_counts = info.get_modification_counts(peptidelist,nterm)
                pep_freq = get_colors.get_position_frequency(protein,peptidelist)
                max_pep_freq = max(pep_freq)

                for key, dict in psm_mod_counts.items():
                        local_max = max(dict.values())
                        if local_max > max_psm:
                            max_psm = local_max



            else:
                pfound = f"""## Showing sample {sample}"""
                peptidelist = protein.get_peptides_by_file_id(sample)
                max_psm = 0
                max_pep_freq = 0

                for s in samples_in_file:
                    samle_peps = protein.get_peptides_by_file_id(s)
                    psm_mod_counts = info.get_modification_counts(samle_peps,nterm)
                    pep_freq = get_colors.get_position_frequency(protein,samle_peps)

                    for key, dict in psm_mod_counts.items():
                        local_max = max(dict.values())
                        if local_max > max_psm:
                            max_psm = local_max
                    if max(pep_freq)>max_pep_freq:
                        max_pep_freq = max(pep_freq)

            mpd =info.get_peptide_modification_dict(peptidelist,selected_modifications,nterm)

            if not peptide_click == None:
                zoomto = [x for x in range(peptide_click[0],peptide_click[1]+1)]
            else:
                zoomto = None
            
            # colorbars = (
            #     html.Label('Peptide abundance colorbar:'),
            #     dcc.Graph(figure = get_colors.create_color_fig1(protein,peptidelist,max_pep_freq)),
            #     html.Label('Residue abundance colorbar:'),
            #     dcc.Graph(figure = get_colors.create_modcolor_fig1(peptidelist,max_psm,selected_modifications))
            # )

            

            ### If the selected modifications are None or the list is empty, set default markdown text, make sure selected modifications is a list, and create the viewer
            if selected_modifications==None or selected_modifications==[] or len(mpd.items())==0:
                text = html.Div([
                            html.H1('Visualizer information'),
                            dcc.Markdown(pfound),
                            #dcc.Markdown(f"""#### Selected modification(s): No modifications have been selected""")
                            dcc.Markdown(f"""#### No modifications have been selected or none of the selected modification types are found""")
                                ])
                selected_modifications=[] # make sure that the selected modifications is an empty list and not a None-value
                
                colorbars = (
                    html.Label('Peptide abundance colorbar:'),
                    dcc.Graph(figure = get_colors.create_color_fig1(protein,peptidelist,max_pep_freq))
                )
                #### If the first MET should be removed, create a viewer that has M1 removed
                if remove_m1:
                    viewer = create_viewer(protein, peptidelist,max_pepfreq_val=max_pep_freq, max_psm=max_psm,nterm=nterm, modifications=selected_modifications, color = viscolor_value, vis_size=vissize_value, labels=label_choice,remove_m1=True,residue_style=resstyle_value,residue_size=ressize_value ,visstyle = visstyle_value, zoomto=zoomto) # No modifications mapped, leader MET removed

                #### Otherwise, create the viewer as is
                else:
                    viewer = create_viewer(protein, peptidelist,max_pepfreq_val=max_pep_freq, max_psm=max_psm,nterm=nterm,modifications=selected_modifications,color = viscolor_value, vis_size=vissize_value, labels=label_choice,remove_m1=False,residue_style=resstyle_value,residue_size=ressize_value,visstyle = visstyle_value, zoomto=zoomto) # No modifications mapped, leader MET present
                
                return view_3d(viewer), text,colorbars # Return the created viewer and the default text

        
            ### If the selected modifications is not empty and not None    
            else:
                merged_mods = info.get_peptide_modification_dict(protein.peptides,selected_modifications,nterm=nterm)    
                merged_mods1 = info.get_peptide_modification_dict(peptidelist,selected_modifications,nterm=nterm) 
                #### If M1 SHOULD be removed
                colorbars = (
                    html.Label('Peptide abundance colorbar:'),
                    dcc.Graph(figure = get_colors.create_color_fig1(protein,peptidelist,max_pep_freq)),
                    html.Label('Residue abundance colorbar:'),
                    dcc.Graph(figure = get_colors.create_modcolor_fig1(peptidelist,max_psm,selected_modifications))
                )
                if remove_m1:

                    new_text = info.get_modinfo_text(selected_modifications,protein,peptidelist,nterm=nterm,remove_m1=True)                                           
                        
                    ##### If only one modification is selected, create the markdown with the new text, DO NOT show frequency plot (since there is only one), and create the viewer without M1
                    if len(selected_modifications)==1:
                        updated_text = html.Div([
                            html.H1('Visualizer information'),
                            dcc.Markdown(pfound),
                            #dcc.Markdown(f"""#### Selected modification(s): {'ation, '.join(selected_modifications)}ation"""),
                            dcc.Markdown(new_text,style={'overflowY':'scroll','height':'500px'})#dash_dangerously_set_inner_html.DangerouslySetInnerHTML(new_text))
                            ])
            
                        viewer = create_viewer(protein, peptidelist,max_pepfreq_val=max_pep_freq, max_psm=max_psm, nterm=nterm,modifications=selected_modifications, color = viscolor_value, vis_size=vissize_value, labels=label_choice,remove_m1 = True,residue_style=resstyle_value,residue_size=ressize_value,visstyle = visstyle_value, zoomto=zoomto) # Only one modification is selected, leader MET is removed
                        
                        return view_3d(viewer), updated_text,colorbars # Return the viewer without M1 and the updated text for the modifications shown
                        
                    ##### If more than one modification is selected, create the markdown with the new text, SHOW the frequency plot, and and create the viewer without M1
                    else:
                        updated_text = html.Div([
                            html.H1('Visualizer information'),
                            dcc.Markdown(pfound),
                            #dcc.Markdown(f"""#### Selected modification(s): {'ation, '.join(selected_modifications)}ation"""),
                            dcc.Markdown(new_text,style={'overflowY':'scroll','height':'500px'})#dash_dangerously_set_inner_html.DangerouslySetInnerHTML(new_text))
                            ])
                        viewer = create_viewer(protein, peptidelist,max_pepfreq_val=max_pep_freq, max_psm=max_psm, nterm=nterm, modifications=selected_modifications,color = viscolor_value, vis_size=vissize_value, labels=label_choice,remove_m1=True,residue_style=resstyle_value,residue_size=ressize_value,visstyle = visstyle_value, zoomto=zoomto) # More than 1 modification is selected, leader MET is removed
                        return view_3d(viewer), updated_text,colorbars # Return the viewer without M1 and the updated text for the modifications shown
                
                #### If M1 should NOT be removed
                else:
                    new_text = info.get_modinfo_text(selected_modifications,protein,peptidelist,nterm=nterm,remove_m1=False)
                    
                    ##### If only one modification is selected, update the text, DO NOT create the frequency plot, create the viewer with M1
                    if len(selected_modifications)==1:
                        updated_text = html.Div([
                            html.H1('Visualizer information'),
                            dcc.Markdown(pfound),
                            #dcc.Markdown(f"""#### Selected modification(s): {'ation, '.join(selected_modifications)}ation"""),
                            dcc.Markdown(new_text,style={'overflowY':'scroll','height':'500px'})#dash_dangerously_set_inner_html.DangerouslySetInnerHTML(new_text))
                        ])


                       
                        viewer = create_viewer(protein,peptidelist,max_pepfreq_val=max_pep_freq, max_psm=max_psm, nterm=nterm, modifications=selected_modifications, color = viscolor_value, vis_size=vissize_value, labels=label_choice,remove_m1 = False,residue_style=resstyle_value,residue_size=ressize_value,visstyle = visstyle_value, zoomto=zoomto) # Only one modification is selected, leader MET is present
                        return view_3d(viewer), updated_text,colorbars # Return the viewer M1 and the updated text for the modifications shown
                        
                    ##### If more than one modification is selected, update the text, create the frequency plot and the viewer with M1
                    else:
                        updated_text = html.Div([
                            html.H1('Visualizer information'),
                            dcc.Markdown(pfound),
                            #dcc.Markdown(f"""#### Selected modification(s): {'ation, '.join(selected_modifications)}ation"""),
                            dcc.Markdown(new_text,style={'overflowY':'scroll','height':'500px'})#dash_dangerously_set_inner_html.DangerouslySetInnerHTML(new_text))
                            ])
                        viewer = create_viewer(protein,peptidelist,max_pepfreq_val=max_pep_freq, max_psm=max_psm, nterm=nterm, modifications=selected_modifications, color = viscolor_value, vis_size=vissize_value, labels=label_choice,remove_m1=False,residue_style=resstyle_value,residue_size=ressize_value,visstyle = visstyle_value, zoomto=zoomto) # More than 1 modification is selected, leader MET is present
                        return view_3d(viewer), updated_text,colorbars # Return the viewer with M1 and the updated text for the modifications shown
            
# @app.callback(
#     Output('met-test-label','children'),
#     Input('remove-met-check','value')
# )
# def met_check(check_val):
#     if check_val:
#         return 'M1 is removed'

@app.callback(
    Output('Click-data','children'),
    Input('pepView','clickData'),
    Input('Zoom-reset','n_clicks'),
    prevent_initial_call = True
)
def get_click_data(clickData,n_clicks):
    """
    Handles the click data from the peptide viewer ('pepView') and zoom reset button ('Zoom-reset')
    to manage zoom levels or other interactive features based on user input.

    Args:
        clickData (dict): Data generated by clicking on the peptide viewer, containing coordinates
                          and other specifics about the point of click.
        n_clicks (int): Number of times the zoom reset button has been clicked.

    Returns:
        list: A list containing the start and end points of the clicked peptide, if 'pepView' was clicked.
              Returns None if the zoom reset button was clicked or if there was no relevant interaction.
    """
    if not ctx.triggered:
        return None

    elif ctx.triggered_id == 'pepView' :
        point_data = clickData['points'][0]
        start = point_data['customdata'][0]
        end = point_data['customdata'][1]
        return [start,end]
    
    elif ctx.triggered_id == 'Zoom-reset':
        return None
    
    

if __name__ == '__main__':
    # port = find_available_port(8051)
    # print(f"Using port: {port}")
    app.run(host='0.0.0.0',port=8051)
