import openpyxl
import pandas as pd
import re
import os

def open_workbook(file_path=''):
    """
    Load a specified Excel workbook and return the 'Proteins' worksheet. The function
    prompts the user to input a file path if none is provided or if the provided file path
    does not exist. Continues to prompt until a valid file path is provided.
    
    Parameters:
    - file_path (str): Path to the Excel file. If empty, the function will prompt for it.
    
    Returns:
    - openpyxl.worksheet.worksheet.Worksheet: Worksheet object for the 'Proteins' sheet.
    
    Example:
    >>> workbook = open_workbook('path/to/file.xlsx')
    >>> print(workbook)
    <Worksheet "Proteins">
    
    Note:
    - This function requires the 'openpyxl' module and 'os' module to be imported.
    - It is assumed that the 'Proteins' sheet exists in the Excel workbook.
    """
    if file_path is None or file_path == '':
        file_path = input("Please provide the filepath below:\n")
    
    while not os.path.exists(file_path):
        file_path = input("Filename not found, please input filename below:\n ")

    workbook = openpyxl.load_workbook(filename=file_path)
    proteins = workbook['Proteins']
    return proteins




def extract_protein_df(protein_workbook):
    """
    Extracts protein data from a specified workbook into a DataFrame, along with indices 
    of proteins found by being marked in the 'Checked' column. 
    
    The function assumes the first row of the worksheet contains column headers and that there is a 'Checked' column indicating which
    rows to include in the resulting DataFrame.

    Parameters:
    - protein_workbook: openpyxl.workbook.workbook.Workbook object containing the protein data.

    Returns:
    - tuple: A DataFrame containing selected protein data based on the 'Checked' column and a list
      of indices for these rows.
    
    Example:
    >>> df, indices = extract_protein_df(workbook)
    >>> print(type(df))
    <class 'pandas.core.frame.DataFrame'>
    >>> print(indices)
    [0, 39, 47, 62, 78]

    Note:
    - Requires the 'pandas' library for DataFrame operations.
    - The function performs data extraction by iterating through rows marked 'Checked' in the 
      workbook and collects these into a DataFrame, while discarding columns entirely filled with 
      'None' values.
    """
    protein_index = []

    data = protein_workbook.values
    cols = next(data)
    df = pd.DataFrame(data, columns=cols)
    for key, item in df['Checked'].items():
        if item:
            protein_index.append(key)

    protein_df = df.iloc[protein_index]
    protein_df = protein_df.dropna(axis=1, how='all')  # Removing columns with all None values
    return protein_df, protein_index





def extract_peptide_df(protein_workbook, accession, proteindf, protein_index):
    """
    Extracts peptide data from a workbook for a specific protein based on its accession number 
    and organizes it into a structured DataFrame. This involves transposing and cleaning the data,
    segmenting it based on protein indices, and applying further processing like extracting specific
    sequences and positions within the proteins.

    Parameters:
    - protein_workbook (Workbook): The workbook from which to extract data, assumed to contain peptide data.
    - accession (str): The accession number of the protein to filter the data for.
    - proteindf (DataFrame): DataFrame containing protein information including accession numbers.
    - protein_index (list of int): List containing the starting index of each protein's data in the DataFrame,
      used to locate and extract relevant peptide data.

    Returns:
    - DataFrame: A DataFrame containing organized peptide information specific to the protein of the
      given accession number. This includes cleaned sequences, positions within the proteins, and modifications.

    Example:
    >>> df = extract_peptide_df(workbook, 'P12345', proteindf, protein_index)
    >>> print(df.head())

    Note:
    - Requires 'pandas' for DataFrame operations and 're' for regular expressions.
    - This function assumes the 'Proteindf' DataFrame includes a column 'Accession' and that the workbook
      data is formatted with peptide information under headers that need specific cleaning and formatting.
    """
    # Changing headers
    df = pd.DataFrame(protein_workbook.values)
    trans = df.T
    trans.pop(0)
    trans.pop(1)
    df = trans.T
    df.pop(0)
    col = df.iloc[0]
    col.name='Peptide index'
    df.columns=col
    df=df[1:]

    # Dividing the dataframe by protein and its peptides
    for i, proteinid in enumerate(protein_index):
        if accession in proteindf['Accession'][proteinid] and not i == len(protein_index)-1:
            pepriderow_start,pepriderow_end = proteinid,protein_index[i+1]-2 #not id+1 and index-2 because dataframe index starts with 1
        elif accession in proteindf['Accession'][proteinid] and i == len(protein_index)-1:
            pepriderow_start,pepriderow_end = proteinid,len(df) #not id+1 and index-2 because dataframe index starts with 1

    newdf = df.iloc[pepriderow_start:pepriderow_end]

    new_df = pd.DataFrame()
    new_df = pd.concat([new_df, newdf], ignore_index=True)

    # Process 'Positions in Proteins', 'Annotated Sequence', and 'Sequence in Protein'
    def extract_positions(text):
        """
        Extracts position information from a text string.

        Parameters:
        - text (str): The text string to extract position information from.

        Returns:
        - tuple or None: A tuple containing the start and end positions if found, otherwise None.

        Example:
        >>> extract_positions('12-34')
        ('12', '34')
        """
        match = re.search(r'(\d+)-(\d+)', text)
        return match.groups() if match else None

    def clean_sequence(sequence, pattern):
        """
        Cleans a sequence string based on a given pattern.

        Parameters:
        - sequence (str): The sequence string to clean.
        - pattern (str): The regular expression pattern to use for cleaning.

        Returns:
        - str or None: The cleaned sequence string if a match is found, otherwise None.

        Example:
        >>> clean_sequence('[AB].[CDE].[FG]', r'\[.*?\]\.(.*?)\.\[.*?\]')
        'CDE'
        """
        match = re.search(pattern, sequence)
        return match.group(1) if match else None

    new_df['Positions in Proteins'] = new_df['Positions in Proteins'].apply(extract_positions)
    new_df['Annotated Sequence'] = new_df['Annotated Sequence'].apply(
        lambda x: clean_sequence(x, r'\[.*?\]\.(.*?)\.\[.*?\]') # lambda x creates an anonymous function that processes each element x in the 'Sequence in Protein' column of new_df, using the clean_sequence function to extract a substring according to the specified regex pattern, and then updates the column with the extracted results.
    )
    new_df['Sequence in Protein'] = new_df['Sequence in Protein'].apply(
        lambda x: clean_sequence(x, r'[^.]+\.(.*?)\.[^.]+') # lambda x creates an anonymous function that processes each element x in the 'Sequence in Protein' column of new_df, using the clean_sequence function to extract a substring according to the specified regex pattern, and then updates the column with the extracted results.
    )
        
    def parse_modifications(mod_string):
        """
        Parses a string representing modifications and organizes them into a list of modification-position pairs.

        Parameters:
        - mod_string (str): The string containing modification information.

        Returns:
        - list: A list of modification-position pairs, where each pair is represented as [modification, position].

        Example:
        >>> parse_modifications('Phospho (STY) (1); Acetyl (K) (6); Methyl (K) (9)')
        [['Phospho', '1'], ['Acetyl', '6'], ['Methyl', '9']]

        Note:
        - Modification information is expected to be provided in the format 'Modification (Position);'.
        - Nested modifications are supported and will be parsed accordingly.
        """
        # Split the string into parts using the semicolon
        parts = mod_string.split(';')
        
        # Create a list to store the modifications
        modifications_list = []

        # Iterate over each part
        for part in parts:
            # Strip whitespace
            clean_part = part.strip()
            if not clean_part:
                continue  # skip empty parts
            
            # Split the part into position and nested modifications
            position_part, nested_mods = clean_part.split('(', 1)
            # Further split nested modifications by closing parenthesis and remove empty strings
            nested_mods = [mod.strip('()') for mod in nested_mods.split(')') if mod]

            # For each modification, append the [modification, position] sublist to modifications_list
            for mod in nested_mods:
                modifications_list.append([mod, position_part.strip()])
    
        return modifications_list

    new_df['Modifications'] = new_df['Modifications'].apply(parse_modifications)

    return new_df

def get_protein(workbook):
    """
    Extracts protein data from a workbook and returns DataFrame, protein indices, and accession numbers.

    Parameters:
    - workbook (Workbook): The workbook containing protein data.

    Returns:
    - tuple: A tuple containing:
        - DataFrame: DataFrame containing protein information.
        - list: List of protein indices.
        - list: List of accession numbers.

    Example:
    >>> protein_df, protein_index, accession_numbers = get_protein(workbook)

    Note:
    - Requires the 'extract_protein_df' function to be defined.
    """
    protein_df,protein_index = extract_protein_df(workbook)

    accession_numbers = []
    for accession in protein_df['Accession']:
        accession_numbers.append(accession)

    return protein_df, protein_index, accession_numbers

def get_samples_in_file(workbook):
    """
    Extracts sample names found in a workbook containing protein data.

    Parameters:
    - workbook (Workbook): The workbook containing protein data.

    Returns:
    - list: List of sample names found in the workbook.

    Example:
    >>> samples_in_file = get_samples_in_file(workbook)

    Note:
    - Requires the 'get_protein' function to be defined.
    """
    protein_df = get_protein(workbook)[0]
    samples_in_file = []
    for column in protein_df:
        if 'Found in Sample' in column:
            # Replace "Found in Sample:" with an empty string and then print
            samples_in_file.append(column.replace("Found in Sample: ", ""))
    return samples_in_file
