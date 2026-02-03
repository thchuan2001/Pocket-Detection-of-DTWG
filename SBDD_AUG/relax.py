from rdkit import Chem
from rdkit.Chem import AllChem
import os
from tqdm import tqdm
import glob
from Bio.PDB import PDBParser, PDBIO
from rdkit.Chem import SDMolSupplier, MolToPDBBlock
from multiprocessing import Pool
import time

SCHRODINGER = '~/schrodinger2024-1'
from Bio import PDB
from Bio.PDB import PDBIO
from rdkit import Chem
from rdkit.Chem import AllChem
import io

def generate_complex_pdb(protein_file, ligand_file, output_file):
    try:
        # Step 1: Read the protein PDB file
        parser = PDB.PDBParser(QUIET=True)
        structure = parser.get_structure('protein', protein_file)
        
        # Merge all chains into a single chain 'A'
        first_model = structure[0]
        new_chain = PDB.Chain.Chain('A')
        for chain in first_model:
            for residue in chain:
                new_chain.add(residue)
        
        # Create a new structure with a single chain
        new_structure = PDB.Structure.Structure('merged_protein')
        new_model = PDB.Model.Model(0)
        new_model.add(new_chain)
        new_structure.add(new_model)
        
        # Step 2: Read the ligand SDF file and convert to PDB format
        ligand_mol = Chem.MolFromMolFile(ligand_file)
        if ligand_mol is None:
            raise ValueError(f"Failed to read ligand from {ligand_file}")
        
        # Convert ligand to PDB format
        ligand_pdb_block = Chem.MolToPDBBlock(ligand_mol)
        
        # Step 3: Read the ligand PDB block and set the chain ID to 'B'
        ligand_structure = parser.get_structure('ligand', io.StringIO(ligand_pdb_block))
        for model in ligand_structure:
            for chain in model:
                chain.id = 'B'
        
        # Step 4: Merge the protein and ligand structures
        for model in ligand_structure:
            for chain in model:
                new_structure[0].add(chain)
        
        # Write the merged structure to the output PDB file
        pdb_io = PDB.PDBIO()
        pdb_io.set_structure(new_structure)
        pdb_io.save(output_file)
    except Exception as e:
        print("Error:", e)
        pass



def relax(dirs,pdbfile,method='local_refine'): 
    start_time = time.time()
    assert method in ['local_refine','mc_refine','minization']
    os.chdir(dirs)
    cmd=f'{SCHRODINGER}/utilities/prepwizard -j prepwizard_{pdbfile.replace(".", "_").replace("/","_")} -watdist 5 -rehtreat -propka_pH 7.4  -HOST localhost:1 -noimpref -TMPLAUNCHDIR -ATTACHED -WAIT {pdbfile} {pdbfile.replace(".pdb","_fixed.maegz")}'

    # print("Running command:", cmd)
    os.system(cmd)
    if method == 'local_refine':
        conf = f'''STRUCT_FILE	{pdbfile.replace(".pdb","_fixed.maegz")}
JOB_TYPE	REFINE
PRIME_TYPE	SITE_OPT
SELECT	asl = fillres within 6.000000 ( res.num 1 AND res.inscode " " AND chain.name B )
LIGAND	asl = res.num 1 AND res.inscode " " AND chain.name B
NPASSES	1
INCRE_PACKING	no
USE_CRYSTAL_SYMMETRY	no
USE_RANDOM_SEED	no
SEED	0
OPLS_VERSION S-OPLS
EXT_DIEL	80.00
USE_MEMBRANE	no
HOST	localhost:1
'''
    elif method == 'mc_refine':
        print('mc_refine is very slow, be patient!')
        conf = f'''STRUCT_FILE	{pdbfile.replace(".pdb","_fixed.maegz")}
JOB_TYPE	REFINE
PRIME_TYPE	MC
SELECT	asl = protein_near_ligand
LIGAND	asl = ligand
NSTEPS	25
PROB_SIDEMC	0.600000
PROB_RIGIDMC	0.300000
PROB_HMC	0.100000
TEMP_SIDEMC	2000.0
TEMP_RIGIDMC	300.0
TEMP_HMC	900.0
FIND_BEST_STRUCTURE	yes
NUM_OUTPUT_STRUCT	1
INCRE_PACKING	no
USE_CRYSTAL_SYMMETRY	no
USE_RANDOM_SEED	no
SEED	0
OPLS_VERSION S-OPLS
EXT_DIEL	80.00
USE_MEMBRANE	no
HOST	localhost:1
'''            
    elif method == 'minization':
        conf = f'''STRUCT_FILE	{pdbfile.replace(".pdb","_fixed.maegz")}
JOB_TYPE	REFINE
PRIME_TYPE	REAL_MIN
SELECT	asl = (fillres within 6.000000 ( res.num 1 AND res.inscode " " AND chain.name B )) or (res.num 1 AND res.inscode " " AND chain.name B)
MINIM_NITER	2
MINIM_NSTEP	65
MINIM_RMSG	0.010000
MINIM_METHOD	au
USE_CRYSTAL_SYMMETRY	no
USE_RANDOM_SEED	no
SEED	0
OPLS_VERSION S-OPLS
EXT_DIEL	80.00
USE_MEMBRANE	no
HOST	localhost:1
'''
    else:
        raise ValueError('method must be one of local_refine,mc_refine,minization')
    with open(pdbfile.replace(".pdb",'_refine.inp'),'w') as f:
        f.write(conf)
    os.system(f'{SCHRODINGER}/prime {pdbfile.replace(".pdb","_refine.inp")} -NJOBS 8 -WAIT')
    #convert to pdb
    os.system(f'{SCHRODINGER}/utilities/structconvert {pdbfile.replace(".pdb","_refine-out.maegz")} {pdbfile.replace(".pdb","_refined.pdb")}')

    # rm temp file
    os.system(f'rm {pdbfile.replace(".pdb","_fixed.maegz")}')
    os.system(f'rm {pdbfile.replace(".pdb","_refine.inp")}')
    os.system(f'rm {pdbfile.replace(".pdb","_refine-out.maegz")}')
    os.system(f'rm {pdbfile.replace(".pdb","_refine.log")}')
    os.system(f'rm prepwizard_{pdbfile.replace(".", "_").replace("/","_")}.log')
    
    time_cost=time.time()-start_time
    return time_cost



# if __name__ == "__main__":

#     dataset=DTWG_DUDEDataset()
#     exp_name="AF2_fpocket_BFN"


#     # generate complex pdb    
#     for item in tqdm(dataset):
#         print(item)
#         if not os.path.exists(item['protein_dir']):
#             continue

#         ligands=glob.glob(os.path.join(item['dir'],exp_name,'*.sdf'))
#         for ligand in ligands:
#             generate_complex_pdb(
#                 protein_file=item['protein_dir'],
#                 ligand_file=ligand,
#                 output_file=ligand.replace('.sdf','_complex.pdb')
#             )

#     # relax pdb 
#     task_list=[]
#     for i,item in tqdm(enumerate(dataset.get_items())):
#         if not os.path.exists(item['protein_dir']):
#             continue
#         complexes=glob.glob(os.path.join(item['dir'],exp_name,'*_complex.pdb'))
#         # sort 
#         complexes.sort()
#         for complex_file in complexes:
#             task_list.append((os.path.dirname(complex_file),os.path.basename(complex_file),'local_refine'))
#     with Pool(101) as p:
#         results = p.starmap(relax, task_list)

def process_single_sdf(protein_file, ligand_file, method='local_refine'):
    """Process a single SDF file: generate complex and relax"""
    try:
        # Generate complex PDB
        output_file = ligand_file.replace('.sdf', '_complex.pdb')
        refined_file = ligand_file.replace('.sdf', '_complex_refined.pdb')
        
        # Skip if already refined
        if os.path.exists(refined_file):
            return ('skip', ligand_file, 0)
        
        # Generate complex if not exists
        if not os.path.exists(output_file):
            generate_complex_pdb(
                protein_file=protein_file,
                ligand_file=ligand_file,
                output_file=output_file
            )
        
        # Relax the complex
        work_dir = os.path.dirname(output_file)
        pdb_filename = os.path.basename(output_file)
        time_cost = relax(work_dir, pdb_filename, method)
        
        return ('success', ligand_file, time_cost)
    except Exception as e:
        return ('error', ligand_file, str(e))


if __name__ == "__main__":

    protein_dir="/home/tanhaichuan/GenPack_for_galaxy/domains"
    bfn_output_dir="/home/tanhaichuan/GenPack_for_galaxy/fpocket_backbone_test_100"

    pocket_list=[f for f in glob.glob(bfn_output_dir+"/*") if os.path.isdir(f)]

    # Collect all SDF files with their corresponding protein files
    task_list = []
    for pocket in tqdm(pocket_list, desc="Collecting tasks"):
        ligand_files = glob.glob(os.path.join(pocket, "*.sdf"))
        protein_file = os.path.join(protein_dir, "_".join(os.path.basename(pocket).split("_")[:3]) + ".pdb")
        
        if not os.path.exists(protein_file):
            print(f"Warning: protein file not found: {protein_file}")
            continue
        
        for ligand_file in ligand_files:
            task_list.append((protein_file, ligand_file, 'local_refine'))
    
    print(f"Total SDF files to process: {len(task_list)}")


    # print the task list
    for task in task_list[:5]:
        print(task)
    print("...")
    
    # Process all SDFs in parallel
    import time as time_module
    start_time = time_module.time()
    
    with Pool(160) as p:
        results = p.starmap(process_single_sdf, task_list)
    
    
    ##########
    # debug
    # task_list = task_list[:1]
    # results = []
    # for task in task_list:
    #     results.append(process_single_sdf(*task))
    ##########
    
    total_time = time_module.time() - start_time
    # Statistics
    success_count = sum(1 for r in results if r[0] == 'success')
    skip_count = sum(1 for r in results if r[0] == 'skip')
    error_count = sum(1 for r in results if r[0] == 'error')
    time_costs = [r[2] for r in results if r[0] == 'success']
    
    print(f"\n=== Processing Results ===")
    print(f"Success: {success_count}")
    print(f"Skipped: {skip_count}")
    print(f"Errors: {error_count}")
    print(f"Total wall time: {total_time:.2f}s ({total_time/60:.2f} minutes)")
    
    if time_costs:
        print(f"\n=== Individual Task Statistics ===")
        print(f"Mean time per task: {sum(time_costs)/len(time_costs):.2f}s")
        print(f"Min time: {min(time_costs):.2f}s")
        print(f"Max time: {max(time_costs):.2f}s")
        print(f"Total processing time (sum): {sum(time_costs):.2f}s ({sum(time_costs)/60:.2f} minutes)")
        print(f"Speedup factor: {sum(time_costs)/total_time:.2f}x")
        
        # Show top 10 slowest tasks
        task_times = [(results[i][1], results[i][2]) for i in range(len(results)) if results[i][0] == 'success']
        task_times.sort(key=lambda x: x[1], reverse=True)
        print(f"\nTop 10 slowest tasks:")
        for i, (ligand, time_cost) in enumerate(task_times[:10]):
            print(f"  [{i+1}] {os.path.basename(ligand)}: {time_cost:.2f}s")
    
    # Print errors
    errors = [(r[1], r[2]) for r in results if r[0] == 'error']
    if errors:
        print(f"\nErrors occurred in {len(errors)} files:")
        for ligand, error in errors[:10]:  # Show first 10 errors
            print(f"  {ligand}: {error}")