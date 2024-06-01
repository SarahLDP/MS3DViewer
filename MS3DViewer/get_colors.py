import matplotlib.pyplot as plt
from matplotlib.colors import to_hex
import numpy as np
from plotly import graph_objects as go

from classes import Protein,Peptide
import info


def get_position_frequency(protein:Protein,peptide_list):
    # Initialize an array to store the frequency of each position
    seq_len = len(protein.master_sequence)
    position_frequency = [0] * (seq_len)


    # Count the frequency of each position
    for peptide in peptide_list:
        peptide:Peptide
        for position in peptide.get_position_range():

            position_frequency[position - 1] += 1
    return position_frequency

def get_peptide_abundance_color(protein:Protein, peptide_list, max_peptide_freq, c_color='cool'):
    freqpos_list = []
    freq_list = []
    pos_freq = get_position_frequency(protein,peptide_list)

    # Adjusting color map to include grey when pos_freq == 0
    min_value = 0
    max_value = max_peptide_freq
    if min_value > 0:
        min_value -= 1
    else:
        min_value -= 0  # no change
    if max_value < 0:
        max_value += 1
    else:
        max_value += 0  # no change

    # Define the color map
    cmap = plt.get_cmap(c_color)

    # Generate colors based on the color map
    colors = {
        i: "#{:02X}{:02X}{:02X}".format(int(255 * color[0]), int(255 * color[1]), int(255 * color[2]))
        for i, color in enumerate(cmap(np.linspace(0, 1, max_value - min_value + 1)))
    }

    # Setting grey color if pos_freq == 0
    if 0 in pos_freq:
        colors[0] = '#808080'

    for i, p in enumerate(pos_freq):
        if i == len(pos_freq) - 1:
            freqpos_list.append(i)
            freq_list.append(p)
        else:
            if p != pos_freq[i + 1]:
                freqpos_list.append(i)
                freq_list.append(p)
    return colors,freqpos_list, freq_list

def set_peptide_abundance_color(colors, viewer,freqpos_list,freq_list,pepstyle = 'cartoon', remove_m1=False,size=0.8):
    selection_list = []

    # if remove_m1:
    #     for i, p in enumerate(freqpos_list):
    #         if i == 0:
    #             pass
    #         elif i == len(freqpos_list):
    #             pass
    #             # selection = f"chain A and :{freqpos_list[i]}-{len(protein.master_sequence)}"
    #         else:
    #             selection = f"chain A and :{freqpos_list[i-1]+1}-{freqpos_list[i]}"
    #             selection_list.append(selection)

    #     print('freqpos_list',freqpos_list)
    #     print('freq_list',freq_list)
    #     print('selection_list',selection_list)

    #     for i, f in enumerate(freq_list):
    #         print('Sel: ',selection_list[i],'     f: ',f,'        color: ',colors[f])
    #         viewer.setStyle({'resi': selection_list[i]}, {pepstyle: {'color': colors[f],'radius':size}})

        
    # else:
    if remove_m1:
        freqpos = []
        for pos in freqpos_list:
            new_pos = pos-1
            freqpos.append(new_pos)
    else:
        freqpos = freqpos_list

    for i, p in enumerate(freqpos):
        if i == 0:
            selection = f"chain A and :{1}-{(freqpos[i]) + 1}"
            selection_list.append(selection)
        if i != 0:
            selection = f"chain A and :{freqpos[i-1] + 2}-{freqpos[i] + 1}"
            selection_list.append(selection)
    

    for i, f in enumerate(freq_list):
        viewer.setStyle({'resi': selection_list[i]}, {pepstyle: {'color': colors[f],'radius':size}})

    return viewer


def get_mod_freq_colors(freq_vals,max_val=None,c_color = 'YlOrRd',isdict=True):
    # Determine the range of values in freq
    if isdict and max_val is None:
        max_value = max(freq_vals.values())
    elif isinstance(freq_vals, (int, float)):
        max_value = freq_vals
    elif max_val is None:
        max_value = max(freq_vals)
    else:
        max_value = max_val

    # Define the color map
    cmap = plt.get_cmap(c_color)
    
    # Generate colors based on the color map
    colors = {i + 1: "#{:02X}{:02X}{:02X}".format(int(255 * color[0]), int(255 * color[1]), int(255 * color[2])) 
              for i, color in enumerate(cmap(np.linspace(0, 1, max_value + 1)))}
    
    return colors


def create_modcolor_fig(colors, modifications=None, title='Color Legend'):
    fig = go.Figure()

    counts = set()
    for position, count in sorted(modifications.items()):
        counts.add(count)
    
    for count in sorted(counts):
        fig.add_trace(go.Scatter(
            x=[1],  # Same x for all to align vertically
            y=[count],
            mode='markers',
            marker=dict(
                color=colors[count],
                size=14,
                line=dict(width=2)
            ),
            name=count,  # Legend name for each marker
            showlegend=True  # Ensure legend entry is shown
        ))

    fig.update_layout(
        title=title,
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(
            showline=False,
            showgrid=False,
            showticklabels=False,
            range = [-1,0]
        ),
        yaxis=dict(
            showline=False,
            showgrid=False,
            showticklabels=False,
            range = [-1,0]
        ),
        showlegend=True,  # Enable legend
        legend=dict(
            x=-1,  # X position of the legend (0 corresponds to the far left)
            y=1.04,  # Y position of the legend (1 corresponds to the top)
            xanchor='left',  # Anchor the legend's x-position from the left
            yanchor='top',    # Anchor the legend's y-position from the top
            title = 'Modification count'
        ),
        margin=dict(l=70, r=70, t=50, b=20),
        width=200,
        height=600,
        title_font = dict(family = 'Georgia',size=18),
        font = dict(family = 'Times New Roman',size=14),
    )
    return fig







def create_color_fig(protein,peptide_list, max_peptide_freq, title = 'Color Legend'):
    colors,freqpos_list, freq_list = get_peptide_abundance_color(protein, peptide_list, max_peptide_freq)

    counts = set()
    for  count in sorted(freq_list):
        counts.add(count)

    fig = go.Figure()
    for count in sorted(counts):
        # Add invisible markers on the plot
        fig.add_trace(go.Scatter(
            x=[1],  # Same x for all to align vertically
            y=[count],
            mode='markers',
            marker=dict(
                color=colors[count],
                size=14,  # Set marker size to 0 to make them invisible on the plot
                line=dict(width=2)
            ),
            name=f'Count: {count}',  # Legend name for each marker
            showlegend=True  # Ensure legend entry is shown
        ))

    fig.update_layout(
        title=title,
        plot_bgcolor='rgba(0,0,0,0)',
        xaxis=dict(
            showline=False,
            showgrid=False,
            showticklabels=False,
            range = [-1,0]
        ),
        yaxis=dict(
            showline=False,
            showgrid=False,
            showticklabels=False,
            range = [-1,0]
        ),
        showlegend=True,  # Enable legend
        legend=dict(
            x=-1,  # X position of the legend (0 corresponds to the far left)
            y=1.04,  # Y position of the legend (1 corresponds to the top)
            xanchor='left',  # Anchor the legend's x-position from the left
            yanchor='top',    # Anchor the legend's y-position from the top
            title = 'Peptide count'
        ),
        margin=dict(l=70, r=70, t=50, b=20),
        width=200,
        height=550
    )

    return fig



def create_color_fig1(protein:Protein,peptide_list, max_peptide_freq, title = 'Color Legend'):
    color_dict,freqpos_list, freq_list = get_peptide_abundance_color(protein, peptide_list,max_peptide_freq)
    if max_peptide_freq==None:
        max_key = max(color_dict.keys())
    else:
        max_key = max_peptide_freq
    color_scale = [(key/max_key, color) for key, color in color_dict.items()]
    heatmap_data = [list(range(max_key+1))]

    fig = go.Figure()

    # Create the heatmap
    fig.add_trace(go.Heatmap(
        z=heatmap_data,
        colorscale=color_scale,
        showscale=False,
    ))

    # Update layout to accommodate the colorbar properly
    fig.update_layout(
        width=400,  # Narrow width, as we only need to show the colorbar
        height=50,  # Adjusted for visibility
        margin=dict(l=10, r=10, t=10, b=10),
        yaxis=dict(
            showticklabels=False,
            showgrid=False,
            showline = False,
        )
    )
 

    return fig


def create_modcolor_fig1(peptide_list, max_psm, modifications=None, title='Color Legend',nterm=False):

    pepmoddict = info.get_peptide_modification_dict(peptide_list,modifications,nterm=nterm)
    # print(len(pepmoddict.items()))
    if len(pepmoddict.items())==0:
        return None
    else:
        if max_psm>max(pepmoddict.values()):
            max_key=max_psm
        else:
            max_key=max(pepmoddict.values())
        # max_key = max(pepmoddict.values())

        cmap = plt.get_cmap('YlOrRd')

        # colors = [cmap(i/max_key) for i in range(max_key)]

        colors = get_mod_freq_colors(freq_vals=pepmoddict,max_val=max_key-1,c_color='YlOrRd')

        hex_colors = ['#FFFFFF']
        hex_colors2 =[to_hex(color) for color in colors.values()]
        hex_colors = hex_colors+hex_colors2

        # print(hex_colors)
        num_colors = len(hex_colors)

        color_scale = [(i / (num_colors - 1), color) for i, color in enumerate(hex_colors)]
            
        fig = go.Figure()

        heatmap_data = [list(range(max_key+2))]

        fig = go.Figure()

        # Create the heatmap
        fig.add_trace(go.Heatmap(
            z=heatmap_data,
            colorscale=color_scale,
            showscale=False,
        ))

        # Update layout to accommodate the colorbar properly
        fig.update_layout(
            width=400,  # Narrow width, as we only need to show the colorbar
            height=50,  # Adjusted for visibility
            margin=dict(l=10, r=10, t=10, b=10),
            yaxis=dict(
                showticklabels=False,
                showgrid=False,
                showline = False,
            )
        )
    

        return fig
