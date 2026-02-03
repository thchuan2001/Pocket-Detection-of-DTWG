import os
from unittest import result
from Bio.PDB import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
import glob
import multiprocessing
from tqdm import tqdm

# def is_illegal_pdb(pdb_file, chain_a_id='A', chain_b_id='B'):
#     parser = PDBParser(QUIET=True)
#     structure = parser.get_structure('structure', pdb_file)

#     # Extract non-hydrogen atoms from chain A and chain B
#     chain_a_atoms = [atom for atom in structure.get_atoms() if atom.get_parent().get_parent().id == chain_a_id and atom.element != 'H']
#     chain_b_atoms = [atom for atom in structure.get_atoms() if atom.get_parent().get_parent().id == chain_b_id and atom.element != 'H']
    
#     # Check for clashes between chain A and chain B atoms
#     neighbor_search = NeighborSearch(chain_a_atoms)
    
#     for atom in chain_b_atoms:
#         neighbors = neighbor_search.search(atom.coord, 2.1)
#         if neighbors:
#             return True
    
#     return False

def nearest_distance(pdb_file, chain_a_id='A', chain_b_id='B'):
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', pdb_file)

        # Extract non-hydrogen atoms from chain A and chain B
        chain_a_atoms = [atom for atom in structure.get_atoms() if atom.get_parent().get_parent().id == chain_a_id and atom.element != 'H']
        chain_b_atoms = [atom for atom in structure.get_atoms() if atom.get_parent().get_parent().id == chain_b_id and atom.element != 'H']
        
        min_distance = float('inf')
        
        for atom in chain_b_atoms:
            for neighbor in chain_a_atoms:
                distance = atom - neighbor
                if distance < min_distance:
                    min_distance = distance
        
        # for atom in chain_a_atoms:
        #     for neighbor in chain_b_atoms:
        #         distance = atom - neighbor
        #         if abs(distance - min_distance) < 1e-6:
        #             print(atom.get_parent().get_parent().id, atom.get_name(), neighbor.get_parent().get_parent().id, neighbor.get_name(), distance)

        return min_distance
    except KeyboardInterrupt:
        raise KeyboardInterrupt
    except Exception as e:
        print (pdb_file)
        print(f'Error: {e}')
        return None

pdb_dir = '/home/tanhaichuan/GenPack_for_galaxy/template_matching_results'
pdb_files = glob.glob(os.path.join(pdb_dir, '*refined.pdb'))

print(f'Found {len(pdb_files)} pdb files')
distance_dist=[]

with multiprocessing.Pool() as pool:
    results=pool.map(nearest_distance, pdb_files)


cnt=0
for pdb_file , result in tqdm(zip(pdb_files, results)):
    if result is None or result < 2.5 or result > 3.5:
        cnt+=1
        os.remove(pdb_file)

print(f'{cnt} illegal pdbs moved')