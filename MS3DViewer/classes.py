
# Standard import
import collections

# Third party imports
from Bio.PDB import PDBParser, Polypeptide

# Custom imports
import url_processing



pdbparser = PDBParser()
ppb = Polypeptide.PPBuilder()


class Protein:
    """
    Represents a protein with methods to retrieve its structural and modification data.
    
    Attributes:
        accession (str): Accession number of the protein.
        total_psms (int): Total number of peptide spectrum matches.
        found_in_sample (dict): Dictionary of the samples as keys and the level of the protein found in the samples as items.
        peptides (list): List of Peptide objects associated with this Protein.
    """
        
    def __init__(self, accession,total_psms):
        self._accession = accession
        self._total_psms = total_psms
        
        self._pdb_file = self.get_pdb_file()
        self._peptides:list[Peptide] = [] # List of peptide objects
        self._master_sequence = self.get_master_sequence()
        

    def __repr__(self) -> str:
        return (f"Protein(accession={self._accession}, total_psms={self._total_psms}")

    def __str__(self) -> str:
        return f"""Protein with accession number: {self.accession}
    * PSM count: {self.total_psms}
    * Protein sequence: {self.master_sequence}
    * PDB file location: {self.pdb_file}
        """
    
    @property
    def accession(self) -> str:
        return self._accession

    @property
    def total_psms(self) -> int:
        return self._total_psms
    
    @property
    def pdb_file(self):
        return self._pdb_file
    
    @property
    def peptides(self):
        return self._peptides
    
    @property
    def master_sequence(self):
        return self._master_sequence
    
    def set_peptides(self,peptidelist):
        self._peptides = peptidelist

    def get_pdb_file(self):
        pdb_file = url_processing.retrieve_fromURL(url_processing.create_url(self.accession))
        return pdb_file
    
    def get_master_sequence(self):
        structure = pdbparser.get_structure(self.accession,self.pdb_file)
        for pp in ppb.build_peptides(structure):
            sequence = pp.get_sequence()
        master_sequence = sequence
        return master_sequence
    
    def get_protein_modification_types(self,nterm:bool):
        modtype_dict = {}

        modification_types = []
        for peptide in self.peptides:
            modified_dict=peptide.get_modified_modification_dict(nterm)
            for key in modified_dict.keys():
                if key not in modification_types:
                    modification_types.append(key)
        modtype_dict[self.accession]=modification_types

        return modtype_dict

    def get_modificationtype_position_frequency(self, mod_of_interest,nterm:bool):
        positions = set()
        for peptide in self.peptides:
            modified_dict=peptide.get_modified_modification_dict(nterm)
            if mod_of_interest in modified_dict.keys():
                for val in modified_dict[mod_of_interest]:
                    positions.add(val)
        positions = list(positions)
        return positions
    
    def get_modification_position_frequency(self, mod_of_interest,nterm:bool):
        positions = []
        for peptide in self.peptides:
            modified_dict=peptide.get_modified_modification_dict(nterm)
            if mod_of_interest in modified_dict.keys():
                for val in modified_dict[mod_of_interest]:
                    positions.append(val)
        
        counter = collections.Counter(positions)
        return counter


    def get_resi_list(self, mod_of_interest,nterm:bool):
        """Returns a list of the positions of the residues with the modification of interest, in the protein"""
        resi_pos_list = []
        resi_mod_pos_freq = self.get_modification_position_frequency(mod_of_interest,nterm)
        for key in resi_mod_pos_freq.keys():
            resi_pos_list.append(key)
        
        return resi_pos_list
    
    def get_peptides_by_file_id(self, fileID):
        peptides = []
        for peptide in self.peptides:
            if peptide.file_id in fileID:
                peptides.append(peptide)
        return peptides

class Peptide:
    """
    Represents a peptide, a component of a protein, with methods to manage its modifications.
    
    Attributes:
        sequence (str): Amino acid sequence of the peptide.
        start_position (int): Start position of the peptide in the protein sequence.
        end_position (int): End position of the peptide in the protein sequence.
        modifications (dict): Dictionary of modifications in the peptide.
        file_id (int): Identifier for the file where this peptide was found.
    """

    def __init__(self, protein:Protein, sequence, start_position, end_position, modifications,file_id):
        self._protein = protein # The associated protein
        self._sequence = sequence
        self._start_position = start_position
        self._end_position = end_position
        self._modifications = modifications
        self._file_id = file_id
        
        self._positions = (self.start_position, self.end_position)
        self._position_range = self.get_position_range()

    def __repr__(self) -> str:
        return (f"Peptide(protein={self._protein}, sequence={self._sequence}, start_position={self._start_position}, end_position={self._end_position}, file_id={self._file_id})")

    def __str__(self) -> str:
        return f"""Peptide assigned to masterprotein: {self._protein._accession}.
    * Peptide sequence: "{self._sequence}".
    * Peptide positions (start, end): {self._positions}.
    * Modifications in peptide: {self._modifications}.
    * File id: {self._file_id}
        """
    
    @property
    def protein(self):
        return self._protein
    
    @property
    def protein_accession(self):
        return self._protein._accession
    
    @property
    def sequence(self):
        return self._sequence

    @property
    def start_position(self):
        return self._start_position

    @property
    def end_position(self):
        return self._end_position

    @property
    def modifications(self):
        return self._modifications

    @property
    def file_id(self):
        return self._file_id


    @property
    def positions(self):
        return self._positions

    @property
    def position_range(self):
        return self._position_range

    def get_position_range(self,remove_M1=False):
        if remove_M1:
            return list(range(self.start_position-1, self.end_position))
        else:
            return list(range(self.start_position, self.end_position + 1))
    
    def get_modified_modification_dict(self,nterm:bool):
        modified_dict = {}
        for modification_type, positions in self.modifications.items():
            if nterm:
                if 'N-Term' in positions:
                    modified_dict[f'{modification_type}|N-Term'] = [self.start_position]
                else:
                    modified_positions = [int(pos[1:]) if pos[1:].isdigit() else 1 for pos in positions]
                    modified_dict[modification_type] = modified_positions
                    modified_positions = [(self.start_position if pos[1:] == 1 else (int(pos[1:]) + self.start_position - 1) if pos[1:] != '' else self.start_position) for pos in positions]
                    modified_dict[modification_type] = modified_positions
            else:
                if 'N-Term' in positions:
                    pass
                else:
                    modified_positions = [int(pos[1:]) if pos[1:].isdigit() else 1 for pos in positions]
                    modified_dict[modification_type] = modified_positions
                    modified_positions = [(self.start_position if pos[1:] == 1 else (int(pos[1:]) + self.start_position - 1) if pos[1:] != '' else self.start_position) for pos in positions]
                    modified_dict[modification_type] = modified_positions
        return modified_dict
