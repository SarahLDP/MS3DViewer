import py3Dmol

from get_colors import get_peptide_abundance_color,set_peptide_abundance_color,get_mod_freq_colors
import info
from classes import Protein
from remove_first_met import remove_first_met

def create_viewer(protein:Protein,peptide_list,nterm, max_pepfreq_val, max_psm, color='By peptide abundance',peptide_coverage_color ='cool',modifications=None,labels='Off',remove_m1=False,residue_style='stick',residue_size = 0.8, visstyle = 'cartoon',vis_size = 0.8,zoomto=None, viewer_height = 800,viewer_width = 1800,mod_freq_color = 'YlOrRd'):
    
    
    if remove_m1:
        file_path = remove_first_met(protein.pdb_file,f'Removed_first_MET-{protein.accession}.pdb')
    else:
        file_path= protein.pdb_file
    

    with open(file_path, 'r') as f:
        pdb_data = f.read()

    viewer = py3Dmol.view(width=viewer_width, height=viewer_height)
    viewer.addModel(pdb_data, 'pdb')
    viewer.setStyle({visstyle: {'colorscheme': {'prop':'resi','min':50,'max':90},'radius':vis_size}})

    if color == 'By peptide abundance':
        if remove_m1:
            peptide_colors,freqpos_list, freq_list = get_peptide_abundance_color(protein, peptide_list,max_pepfreq_val,peptide_coverage_color)            
            set_peptide_abundance_color(peptide_colors,viewer,freqpos_list,freq_list, pepstyle = visstyle, remove_m1=True,size=vis_size)

        else:
            peptide_colors,freqpos_list, freq_list = get_peptide_abundance_color(protein,peptide_list,max_pepfreq_val,peptide_coverage_color)
            set_peptide_abundance_color(peptide_colors,viewer,freqpos_list,freq_list, pepstyle = visstyle, remove_m1=False,size=vis_size)

    elif color == 'Rainbow':
        viewer.setStyle({visstyle: {'colorscheme': {'prop':'resi','gradient': 'roygb','min':50,'max':90},'radius':vis_size}})
    else:
        viewer.setStyle({visstyle: {'colorscheme': 'gray','radius':vis_size}})
 
    if not modifications==[] and not modifications==None:
        if remove_m1:
            
            mpd =info.get_peptide_modification_dict(peptide_list,modifications,nterm)
            if len(mpd.items())==0:
                pass
            else:
                if max_psm>max(mpd.values()):
                    max_val=max_psm
                else:
                    max_val=max(mpd.values())

                colors = get_mod_freq_colors(freq_vals=mpd, max_val=max_val,c_color=mod_freq_color)

                for mod_position,mod_count in mpd.items():
                    viewer.addStyle({"chain": 'A', "resi": mod_position-1},
                                            {residue_style: {"color": colors[mod_count], "radius": residue_size}})
                    if labels=='On':
                        viewer.addResLabels({"chain": 'A',"resi": mod_position-1},{"backgroundColor": "lightgray","fontColor": "black","backgroundOpacity": 0.1})

        else:
            mpd =info.get_peptide_modification_dict(peptide_list,modifications,nterm)

            if len(mpd.items())==0:
                pass
            else:
            
                if max_psm>max(mpd.values()):
                    max_val=max_psm
                else:
                    max_val=max(mpd.values())
                
                colors = get_mod_freq_colors(freq_vals=mpd,max_val=max_val,c_color=mod_freq_color)
                
                for mod_position,mod_count in mpd.items():
                    viewer.addStyle({"chain": 'A', "resi": mod_position},
                                            {residue_style: {"color": colors[mod_count], "radius": residue_size}})
                    if labels=='On':
                        viewer.addResLabels({"chain": 'A',"resi": mod_position},{"backgroundColor": "lightgray","fontColor": "black","backgroundOpacity": 0.1})
    
    #viewer.addSurface(py3Dmol.SAS,{'opacity':0.5,'color':'blue'})  
    if zoomto == None:
        viewer.zoomTo()
    else:
        viewer.zoomTo({'resi':zoomto})
        viewer.addSurface(py3Dmol.SAS,{'opacity':0.5,'color':'white'},{'resi':zoomto})

    return viewer
