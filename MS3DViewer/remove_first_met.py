from Bio.PDB import PDBParser,PDBIO
import os

parser=PDBParser()

def remove_first_met(pdb_file,filename = None):

    if os.path.exists(filename):
        return filename
    
    else:
        structure = parser.get_structure('pdb_structure', pdb_file)

        # Iterate over chains in the structure
        for chain in structure.get_chains():
            residues = list(chain.get_residues())

            # Check if the first residue is MET
            if residues[0].get_resname() == "MET":
                # Remove the first residue
                chain.detach_child(residues[0].id)

                # # Update residue indices
                # for i, residue in enumerate(chain.get_residues(), start=1):
                #     residue.id = (' ', i, ' ')

            for i, residue in enumerate(chain.get_residues(), start=1):
                # Handling of the residue ID which includes hetero-flag, sequence number, and insertion code
                hetfield, resseq, icode = residue.id
                residue.id = (hetfield, i, icode)


        # Save the modified structure
        if filename == None:
            filename = "Removed_First_MET"
            io = PDBIO()
            io.set_structure(structure)
            io.save(filename)
        else:
            io = PDBIO()
            io.set_structure(structure)
            io.save(filename)
        return filename
