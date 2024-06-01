import os
import requests


def create_url(accession_number,type='pdb',source='alphafold'):
    """
    Creates the url to the accession number of interest with the filetype of interest
    Type = xml(.xml file from uniprot), fasta(.fasta file from uniprot), alphafold((.pdb file from alphafold with predicted structure))
    """
    uniprot_xml_url = f"https://www.uniprot.org/uniprotkb/{accession_number}.xml"
    uniprot_fasta_url = f"https://www.uniprot.org/uniprotkb/{accession_number}.fasta"
    alphafold_pdb_url = f"https://alphafold.ebi.ac.uk/files/AF-{accession_number}-F1-model_v4.pdb"
    alphafold_cif_url = f'https://alphafold.ebi.ac.uk/files/AF-{accession_number}-F1-model_v4.cif'

    if type == 'xml':
        if source == 'uniprot':
            url = uniprot_xml_url
    elif type == 'fasta':
        if source == 'uniprot':
            url = uniprot_fasta_url
    elif type == 'pdb':
        if source=='alphafold':
            url = alphafold_pdb_url
    elif type == 'mmcif' or type == 'cif':
        if source == 'alphafold':
            url = alphafold_cif_url
    else:
        print('Not valid type')

    return url

def retrieve_fromURL(url, folder= 'Current Folder'):
    """
    Retrieves the files from the url of interest to the folder of choice, and returns the filename.
    """
    if folder == 'Current Folder':
    # Specify folder to save files in
        model_folder = os.getcwd() #os.getcwd() chooses current folder
    # Make sure the folder exists
        os.makedirs(model_folder, exist_ok=True)
        folder = model_folder
    
    response = requests.get(url)
    if response.status_code == 200:
        data = response.text
        # Extract the filename from the URL
        filename = os.path.join(folder, os.path.basename(url))
        # Check if the file already exists
        if os.path.exists(filename):
            pass
        else:
            # Save the data to the specified folder
            with open(filename, 'w') as file:
                file.write(data)
            print(f"Downloaded {filename}")
    else:
        print("URL not found")
    return filename