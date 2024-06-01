# MS3DViewer
MS3DViewer is a comprehensive bioinformatics toolkit designed to facilitate the analysis and 3D visualization of protein structures and their related peptide data. Built with Python, it provides a suite of tools that enable researchers to handle protein data efficiently, from parsing and extraction to advanced visual representation in 3D.

## Features
Data Parsing: Extract protein and peptide data from various formats, including Excel and external databases like UniProt and AlphaFold.
Protein Modifications: Analyze and process different types of protein modifications to better understand functional implications.
3D Visualization: Utilize py3Dmol to visualize protein structures in three dimensions with options to highlight modifications and peptide coverage.
Customizability: Scripts support customization of visual outputs to suit specific research needs, emphasizing flexibility in data presentation.

## Components
classes.py: Defines classes for managing protein and peptide data.
get_colors.py: Tools for assigning colors based on peptide data.
info.py, parse_file.py: Modules for extracting and processing protein data from files.
MS3Dviewer.py: Main script for 3D visualization of proteins in web applications.
peptide_atlas.py: Generates detailed visualizations of peptide mappings on proteins.
remove_first_met.py: Handles the removal of the initial methionine from protein sequences for structural studies.
symbol_assignation.py: Assigns symbols and colors to different protein modifications for visualization.
url_processing.py: Retrieves protein data files from specified URLs.
viewer.py: Configures and displays 3D visualizations of protein structures.

## Getting Started
**1. Clone the Repository:**

git clone https://github.com/yourusername/MS3DViewer.git

2. Install Dependencies:
Make sure Python is installed on your system and install the required libraries:

pip install -r requirements.txt

3. Run the MS3DViewer:
Navigate to the directory containing the scripts and run the desired Python scripts as follows:

python MS3DViewer.py

## Contributing
Contributions to MS3DViewer are welcome! Please fork the repository and submit a pull request with your changes. For major changes, please open an issue first to discuss what you would like to change.

## License
This project is licensed under the MIT License
