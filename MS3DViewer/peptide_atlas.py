
import pandas as pd
import plotly.express as px
import ast

from symbol_assignation import symbol_assignation
import get_colors
from classes import Peptide

def peptide_atlas(peptides:list[Peptide],nterm,gap:int=1,range_max = None,remove_m1=False):
    if remove_m1:
        data = {
            "Sequence": [peptide.sequence for peptide in peptides],
            "Proteins": [peptide.protein_accession for peptide in peptides],
            "Start position": [peptide.start_position-1 for peptide in peptides],
            "End position":[peptide.end_position-1 for peptide in peptides],
            "Position range":[peptide.get_position_range(remove_M1=True) for peptide in peptides],
            "Modifications":[peptide.modifications for peptide in peptides],
            "Modification position": [peptide.get_modified_modification_dict(nterm) for peptide in peptides]
            }
        sequence = peptides[0].protein.master_sequence[1:]
    else:
        data = {
            "Sequence": [peptide.sequence for peptide in peptides],
            "Proteins": [peptide.protein_accession for peptide in peptides],
            "Start position": [peptide.start_position for peptide in peptides],
            "End position":[peptide.end_position for peptide in peptides],
            "Position range":[peptide.get_position_range() for peptide in peptides],
            "Modifications":[peptide.modifications for peptide in peptides],
            "Modification position": [peptide.get_modified_modification_dict(nterm) for peptide in peptides]
            }
        sequence = peptides[0].protein.master_sequence

    pg = pd.DataFrame(data).sort_values(['Start position', 'End position'])

    grouppg = pg
    # Convert columns containing lists or dictionaries to strings
    grouppg['Position range'] = grouppg['Position range'].astype(str)
    grouppg['Modifications'] = grouppg['Modifications'].astype(str)
    grouppg['Modification position'] = grouppg['Modification position'].astype(str)

    # Define a function to safely convert string back to original type
    def convert_to_original(data):
        try:
            return ast.literal_eval(data)
        except (ValueError, SyntaxError):
            return data  # Return as is if conversion is not possible


    # Group by all columns and count the occurrences
    grouped = grouppg.groupby(list(grouppg.columns)).size().reset_index(name='Count')

    grouped['Position range'] = grouped['Position range'].apply(convert_to_original)
    grouped['Modifications'] = grouped['Modifications'].apply(convert_to_original)
    grouped['Modification position'] = grouped['Modification position'].apply(convert_to_original)

    grouped = grouped.sort_values(['Start position', 'End position'])

    grouped.to_excel('peptide_atlas.xlsx',index=False)

    
    trace_ends, traces = [-gap], []
    for s,e in zip(grouped["Start position"],grouped["End position"]):
        traced = False
        for t in range(len(trace_ends)):
            if s >= trace_ends[t]+gap:
                traces.append(t)
                trace_ends[t] = e
                traced = True
                break
        if not traced:
            traces.append(len(trace_ends))
            trace_ends.append(e)
    grouped.insert(0,"Trace",traces)
    grouped.insert(0,"Index",grouped.index)


    x_v = [x for x in grouped['Position range']]
    y_v = [y for y in grouped["Trace"]]
    ps = [s for s in grouped['Start position']]
    pe = [e for e in grouped["End position"]]
    c = [c for c in grouped['Count']]

    
    x_values = []
    y_values = []
    pep_index = []
    pep_start = []
    pep_end = []
    counts = []


    
    
    for idx,x_list in enumerate(x_v):
        for x in x_list:
            x_values.append(x)
            y_values.append(y_v[idx])
            pep_start.append(ps[idx])
            pep_end.append(pe[idx])
            pep_index.append(idx)
            counts.append(c[idx])


    # Create DataFrame
    df = pd.DataFrame({'Residue': x_values, 
                       'y_values': y_values,
                       'Start position':pep_start,
                       'End position':pep_end,
                       'pep_index':pep_index,
                       'PSM Count':counts

                      })

    colors = get_colors.get_mod_freq_colors(freq_vals=counts,c_color='viridis_r',isdict=False)

    set_colors = [colors[count] for count in grouped['Count']]

    pa = px.line(df,
                 x="Residue", 
                 y="y_values", 
                 template="ggplot2", 
                 color="pep_index",
                 height=100+30*(max(y_values)+3),
                 range_x=[0,max(x_values) if range_max is None else range_max],
                 range_y=[-0.5,(max(y_values)+0.5)],
                 color_discrete_sequence=set_colors,
                 title=f"Peptides identified for {grouped.Proteins.iloc[0]}",
                 hover_data = {'Residue':True,'y_values':False,'Start position':True,'End position':True,'pep_index':False,'PSM Count':True})\
    .update_traces(line_width=5,showlegend=False).update_yaxes(visible=False).update_xaxes(showgrid=True,title_text = '')

    
    ticktext = []
    for letter in sequence:
        ticktext.append(letter)
    tickvals = list(range(1,len(sequence)+1))
    pa.update_xaxes(ticktext=ticktext,tickvals=tickvals,side='top')

    pa.update_layout(
        title_font=dict(family='Georgia',size = 18),
    )
    return pa, pg, grouped




def get_modification_scatter_values(pg,mod_of_interest):
    x_values = []
    y_values = []
    start = []
    end = []

    for index,row in pg.iterrows():
        modifications = row["Modification position"]
        for modification,positions in modifications.items():
            if modification == mod_of_interest:
                for position in positions:
                    x_values.append(position)
                    trace_value = pg.at[index,'Trace']
                    y_values.append(trace_value)
                    peptide_start = pg.at[index,'Start position']
                    peptide_end = pg.at[index,'End position']
                    start.append(peptide_start)
                    end.append(peptide_end)
    return x_values,y_values,start,end


def add_modification_legend(pa, pg, protein,nterm, remove_m1=False):
    for modis in protein.get_protein_modification_types(nterm).values():
        for modi in modis:
            if modi=='Phospho':
                modification='Phosphorylation'
            else:
                modification=f'{modi}ation'
            #symbol depending on modification
            symbol,color,linecolor,size = symbol_assignation(modi)
                
            x_val, y_val, start, end = get_modification_scatter_values(pg, modi)
            custom_data = [(s, e) for s, e in zip(start, end)]
            if remove_m1:
                new_x = [x-1 for x in x_val]
                pa.add_scatter(
                    showlegend=True, 
                    x=new_x, 
                    y=y_val,
                    customdata = custom_data, 
                    mode='markers', 
                    legendgroup=modi,
                    name=modi, 
                    marker=dict(
                        color=color,
                        size=size,
                        symbol=symbol,
                        line=dict(
                            color=linecolor,
                            width = 2
                        )
                    ),
                    hovertemplate=f"<b>{modification}</b><br>" +"Position in protein: %{x}%{x:,.0f}<br>"+"<extra></extra>",
                )
                #color_index = (color_index + 1) % len(colors)

            else:
                pa.add_scatter(
                    showlegend=True, 
                    x=x_val, 
                    y=y_val,
                    customdata = custom_data,
                    mode='markers', 
                    legendgroup=modi, 
                    name=modi, 
                    marker=dict(
                        color=color,
                        size=size,
                        symbol=symbol,
                        line=dict(
                            color=linecolor,
                            width = 2
                        )
                    ),
                    hovertemplate=f"<b>{modification}</b><br>" +"Position in protein: %{x:,.0f}<br>"+"<extra></extra>",
                )
                #color_index = (color_index + 1) % len(colors)
    pa.update_layout(legend_title_text = 'Modifications',legend_title_font=dict(family='Georgia',size = 18),legend_font=dict(family='Georgia',size = 14))
    return pa
