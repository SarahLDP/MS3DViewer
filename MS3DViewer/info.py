from collections import OrderedDict, Counter, defaultdict

# Custom module imports
import parse_file
from classes import Protein,Peptide

def get_protein(workbook):
    protein_df,protein_index = parse_file.extract_protein_df(workbook)

    accession_numbers = []
    for accession in protein_df['Accession']:
        accession_numbers.append(accession)

    return protein_df, protein_index, accession_numbers

def get_samples_in_file(workbook):
    protein_df = get_protein(workbook)[0]
    samples_in_file = []
    for column in protein_df:
        if 'Found in Sample' in column:
            # Replace "Found in Sample:" with an empty string and then print
            samples_in_file.append(column.replace("Found in Sample: ", ""))
    return samples_in_file


def merge_dicts(dicts):
    merged_dict = OrderedDict()
    for d in sorted(dicts):
        for key, value in sorted(d.items()):
            if key in merged_dict:
                # If the value is numeric, add them together
                if isinstance(value, (int, float)):
                    merged_dict[key] += value
                else:
                    # If the value is not numeric, just update
                    merged_dict[key] = value
            else:
                merged_dict[key] = value
    return merged_dict

def get_merged_modifications(nterm, modifications:list,protein:Protein,remove_m1=False):
    mod_pos_list = []
    for mod in modifications:
        residues_of_interest = protein.get_resi_list(nterm,mod)
        if residues_of_interest == []:
            pass
        else:
            if remove_m1:
                mod_pos = protein.get_modificationtype_position_frequency(mod)
                new_mod_pos = [m-1 for m in mod_pos]
                mod_pos_list.append(new_mod_pos)
            else:
                mod_pos_list.append(protein.get_modificationtype_position_frequency(mod))

    mod_freq_list = [Counter(l) for l in mod_pos_list]
    merged_mods = merge_dicts(mod_freq_list)
    return merged_mods

def get_peptidemod_counts(peptide_list,nterm):
    """counts the amount of times a modification is found at a specific position"""
    moddictlist = []
    for peptide in peptide_list:
        peptide:Peptide
        moddict = peptide.get_modified_modification_dict(nterm)
        moddictlist.append(moddict)
    
    # Dictionary to count modifications at specific positions
    modification_counts = defaultdict(lambda: defaultdict(int))

    # Count each modification at each position
    for mod_dict in moddictlist:
        for modification, positions in mod_dict.items():
            for position in positions:
                modification_counts[modification][position] += 1
    
    return modification_counts

def get_modification_counts(peptide_list,nterm):
    
    moddictlist = []
    for peptide in peptide_list:
        peptide:Peptide
        moddict = peptide.get_modified_modification_dict(nterm)
        moddictlist.append(moddict)
 
    # Dictionary to count modifications at specific positions
    modification_counts = defaultdict(lambda: defaultdict(int))

    # Count each modification at each position
    for mod_dict in moddictlist:
        for modification, positions in mod_dict.items():
            for position in positions:
                modification_counts[modification][position] += 1

    return modification_counts

def get_modinfo_text(selected_modifications,protein:Protein,peptide_list,nterm,remove_m1=False):
    new_text = '' # Make an empty string for the new informaiton text

    # moddictlist = []
    # for peptide in peptide_list:
    #     peptide:Peptide
    #     moddict = peptide.get_modified_modification_dict(nterm)
    #     moddictlist.append(moddict)
 
    # # Dictionary to count modifications at specific positions
    # modification_counts = defaultdict(lambda: defaultdict(int))

    # # Count each modification at each position
    # for mod_dict in moddictlist:
    #     for modification, positions in mod_dict.items():
    #         for position in positions:
    #             modification_counts[modification][position] += 1
    modification_counts=get_modification_counts(peptide_list,nterm)

    for modification, positions in modification_counts.items():
        for s_mod in selected_modifications:
            if s_mod == modification:
                if modification == 'Phospho':
                    new_text+= f"""### Modification: {modification}rylation.\n"""
                else:
                    new_text+= f"""### Modification: {modification}ation.\n"""
                
                sorted_positions = sorted(positions)
                for position in sorted_positions:
                   
                    if remove_m1:
                        new_text+= f"""- Found at : {protein.master_sequence[position-1]}{position-1} in {modification_counts[modification][position]} PSMs.\n"""
                    else:
                        new_text+= f"""- Found at : {protein.master_sequence[position-1]}{position} in {modification_counts[modification][position]} PSMs.\n"""

    return new_text

def get_peptide_modification_dict(peptide_list,modifications,nterm):
    modcounts = get_peptidemod_counts(peptide_list,nterm)
    list_of_moddicts = []
            
    for modification, position in modcounts.items():
        for mod in modifications:
            if mod == modification:
                list_of_moddicts.append(position)
    
    mod_pos_dict = defaultdict(int)
    
    for d in list_of_moddicts:
        for key,value in d.items():
            mod_pos_dict[key] += value
    
    return mod_pos_dict
