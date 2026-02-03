import os
import sys
from tqdm import tqdm
from rdkit import Chem
import numpy as np
from multiprocessing import Pool
from Bio.PDB import PDBParser, PDBIO
import glob

class PocketExtractor():
    
    def __init__(self):
        pass

    def run(self, tasks):
        # task: {protein:"", ligand:"", threshold:"",output:""}
        # Use pool.map to process tasks in parallel
        task_list=[(x['protein'],x['ligand'],x['threshold'],x['output']) for x in tasks]
        with Pool(128) as p:
            results = p.starmap(self._extract_single, task_list)
    
    def _extract_single(self,protein_file,ligand_file,threshold,output_file):

        if ligand_file is not None:
            # read ligand
            try:
                if ligand_file.endswith('.mol2'):
                    ligand = Chem.MolFromMol2File(ligand_file,sanitize=False)
                elif ligand_file.endswith('.pdb'):
                    ligand = Chem.MolFromPDBFile(ligand_file,sanitize=False)
                elif ligand_file.endswith('.sdf'):
                    ligand = Chem.MolFromMolFile(ligand_file,sanitize=False)
                else:
                    raise NotImplementedError
                assert ligand is not None
            except:
                print(f"Failed to read ligand {ligand_file}")
                return
                
            conf = ligand.GetConformer()
            ligand_coords = conf.GetPositions()
            
            # read protein
            try:
                protein = PDBParser(QUIET=True).get_structure("protein",protein_file)[0]
                assert protein is not None
            except:
                print(f"Failed to read protein {protein_file}")
                return
        else:
            # read protein (complex with chain A and B)
            try:
                protein = PDBParser(QUIET=True).get_structure("protein",protein_file)[0]
                assert protein is not None
            except:
                print(f"Failed to read protein {protein_file}")
                return

            # Extract ligand coordinates from chain B
            if 'B' not in protein:
                print(f"No chain B found in {protein_file}")
                return
                
            ligand_chain = protein['B']
            ligand_coords = []
            for residue in ligand_chain:
                for atom in residue:
                    ligand_coords.append(atom.coord)

            # Remove chain B, keep only chain A
            for chain in list(protein):
                if chain.id != 'A':
                    protein.detach_child(chain.id)

        # extract pocket (residues within threshold distance of ligand)
        for chain in protein:
            remove_residue_ids=[]
            for residue in chain:
                f=1
                for atom in residue:
                    protein_atom_coords = np.array(atom.coord)
                    for ligand_coord in ligand_coords:
                        if np.linalg.norm(protein_atom_coords - ligand_coord) < threshold:
                            f=0
                            break
                if f:
                    remove_residue_ids.append(residue.id)
            for residue_id in remove_residue_ids:
                chain.detach_child(residue_id)
            
        self.pocket = protein
        # remove all atoms with element X or H
        self._remove_atom_element_X_H()

        # save pocket
        io = PDBIO()
        io.set_structure(protein)
        io.save(output_file)
    
    def _remove_atom_element_X_H(self):
        for chain in self.pocket:
            for residue in chain:
                remove_atom_ids = []
                for atom in residue:
                    if atom.element == 'X':
                        remove_atom_ids.append(atom.id)
                    if atom.element == 'H':
                        remove_atom_ids.append(atom.id)
                for atom_id in remove_atom_ids:
                    residue.detach_child(atom_id)


if __name__ == "__main__":

    input_dir = "/home/tanhaichuan/GenPack_for_galaxy/template_matching_results"
    
    # Find all *_refined.pdb files
    print("Searching for refined complex files...")
    complexes = glob.glob(os.path.join(input_dir, '*_refined.pdb'))
    print(f"Found {len(complexes)} refined complex files")
    
    # extract pocket
    tasks = []
    for complex_file in tqdm(complexes, desc="Preparing tasks"):
        # Skip if pocket file already exists
        output_file = complex_file.replace('_refined.pdb', '_refined_pocket6A.pdb')
        if os.path.exists(output_file):
            continue
            
        tasks.append({
            'protein': complex_file,
            'ligand': None,  # Ligand is in chain B of the complex
            'threshold': 6,
            'output': output_file
        })

    print(f"Processing {len(tasks)} tasks...")
    if len(tasks) > 0:
        pocket_extractor = PocketExtractor()
        pocket_extractor.run(tasks)
        print(f"Completed! Generated {len(tasks)} pocket files")
    else:
        print("No new pockets to extract (all already exist)")
