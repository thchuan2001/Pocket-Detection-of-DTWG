import os
import sys
import glob
from tqdm import tqdm
from Bio.PDB import PDBParser
from multiprocessing import Pool
import warnings
from Bio.PDB.StructureBuilder import PDBConstructionWarning

# Import functions from utils.py
sys.path.insert(0, '/home/tanhaichuan/GenPack_for_galaxy/Pocket-Detection-of-DTWG/utils')
from utils import pdb2dict, write_lmdb

warnings.filterwarnings(action='ignore', category=PDBConstructionWarning)


def process_single_pocket(pocket_file):
    """Process a single pocket file - to be used in parallel"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('pocket', pocket_file)
        pocket_name = pocket_file
        dic = pdb2dict(structure, pocket_name)
        return dic
    except Exception as e:
        print(f"Failed to process {pocket_file}: {e}")
        return None


def main():
    # template_matching_results directory
    template_matching_dir = "/home/tanhaichuan/GenPack_for_galaxy/template_matching_results"
    
    # Output LMDB file
    output_lmdb = os.path.join("/home/tanhaichuan/GenPack_for_galaxy/template_matching_pockets.lmdb")
    
    # Step 1: Glob all pocket files first
    print("Searching for all pocket files...")
    all_pocket_files = glob.glob(os.path.join(template_matching_dir, "*_pocket6A.pdb"))
    
    print(f"Found {len(all_pocket_files)} pocket files")
    
    if len(all_pocket_files) == 0:
        print("No pocket files found!")
        return
    
    # Step 2: Process all pocket files in parallel using multiprocessing
    print(f"Processing {len(all_pocket_files)} pocket files in parallel...")
    num_processes = 64  # Adjust based on your system
    
    with Pool(processes=num_processes) as pool:
        results = list(tqdm(
            pool.imap(process_single_pocket, all_pocket_files),
            total=len(all_pocket_files),
            desc="Processing pockets"
        ))
    
    # Filter out None values (failed processing)
    dics = [r for r in results if r is not None]
    
    print(f"Successfully processed {len(dics)} pocket files")
    failed_count = len(all_pocket_files) - len(dics)
    if failed_count > 0:
        print(f"Failed to process {failed_count} files")
    
    # Step 3: Write to LMDB
    if dics:
        print(f"Writing {len(dics)} entries to LMDB...")
        write_lmdb(dics, output_lmdb)
        print(f"Successfully wrote to {output_lmdb}")
    else:
        print("No pocket files were successfully processed!")


if __name__ == "__main__":
    main()
