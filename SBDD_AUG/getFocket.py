import os
from Bio.PDB import PDBParser,Chain,Model,Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import is_aa
from Bio.PDB.Residue import DisorderedResidue,Residue
from Bio.PDB.Atom import DisorderedAtom
import warnings
from Bio.PDB.StructureBuilder import PDBConstructionWarning
import numpy as np
from rdkit.Chem import Descriptors
import copy
import numpy as np
#from biopandas.mol2 import PandasMol2
from rdkit import Chem
import rdkit.Chem.AllChem as AllChem
from tqdm import tqdm
from io import StringIO
import multiprocessing as mp
import freesasa
from copy import deepcopy

warnings.filterwarnings(
    action='ignore',
    category=PDBConstructionWarning)

def get_binding_pockets(original_pdb,lig_coord,thres=10,rm_sidchain=True):

    chain = original_pdb[0]['A'] #only deal with A chain
    tmp_chain = Chain.Chain('A') 
    resid = set()
    for res in chain:
        res_coord = np.array([i.get_coord() for i in res.get_atoms() if i.element != 'H'])
        dist = np.linalg.norm(res_coord[:,None,:]-lig_coord[None,:,:],axis=-1).min()
        if dist<=thres:
            tmp_chain.add(res.copy())
            resid.add(str(res.id[1]))
    tmp_structure = Structure.Structure(original_pdb.id)
    tmp_model = Model.Model(0)
    tmp_structure.add(tmp_model) 
    tmp_model.add(tmp_chain)
    if rm_sidchain:
        for res in tmp_chain:
            remove_atom_ids = []
            for atom in res.get_atoms():
                if not atom.name in ['C','CA','CB','N','O']:
                    remove_atom_ids.append(atom.id) #should not iter over a modified collection
            for atom_id in remove_atom_ids:
                res.detach_child(atom_id)
    return tmp_structure,resid

def pqr_parser(filename,return_score=False):
    with open(filename,'r') as f:
        data = f.readlines()
    coord = []
    for l in data:
        if "Drug Score" in l:
            score = float(l.split()[-1])
        elif l[:4] == 'ATOM':
            coord.append(
                [float(l[30:38]),float(l[38:46]),float(l[46:54])]
            )
    coord = np.array(coord)
    if return_score:
        return coord,score
    else:
        return coord

def extract_pocket(pdb,output_dir,dist_thres=10,rm_sidchain=True,score_thres=0.001,size_thres=35):
    p = PDBParser()
    pdb_struct = p.get_structure('protein',pdb)
    pockets_dir = pdb.replace('.pdb','_out/pockets')
    os.chdir(os.path.dirname(pdb))
    if not os.path.isdir(pockets_dir):
        os.system(f"fpocket -f {pdb} > /dev/null 2>&1")
    pocket_list = sorted(os.listdir(pockets_dir),key=lambda x:int(x.split('_')[0][6:]))
    res_id_list = []
    for f in pocket_list:
        if f[-4:] == '.pqr':
            pocket_id = f.split('_')[0][6:]
            lig_coord,score = pqr_parser(os.path.join(pockets_dir,f),return_score=True)
            if score>=score_thres:
                pocket_name = os.path.basename(pdb).replace(".pdb",f"_{pocket_id}.pdb")
                pocket_structure,resid = get_binding_pockets(pdb_struct,lig_coord,dist_thres,rm_sidchain)
                if len(resid)>=size_thres:
                    io = PDBIO()
                    io.set_structure(pocket_structure)
                    io.save(os.path.join(output_dir,pocket_name))
                    res_id_list.append(f"{pocket_name}\t{score}\t{len(resid)}\t"+','.join(list(resid))+'\n')
    return res_id_list

def process_single_pdb(pdbid, pdb_dir, output_dir):
    """Process a single PDB file and return the results"""
    try:
        results = extract_pocket(
            pdb=os.path.join(pdb_dir, pdbid),
            output_dir=output_dir,
        )
        return ('success', pdbid, results)
    except Exception as e:
        return ('error', pdbid, str(e))

if __name__ == '__main__':
    import random
    import multiprocessing
    from tqdm import tqdm
    from functools import partial

    pdb_dir = '/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads/domains'
    output_dir = '/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads/domains_pocket_backbone'

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    N = 64  # 使用64个并行进程
    log = open('fpocket_detection.log','w')

    pdb_list = [i for i in os.listdir(pdb_dir) if i[-4:] == '.pdb']
    random.seed(42)
    random.shuffle(pdb_list)
    
    # 使用多进程池处理
    process_func = partial(process_single_pdb, pdb_dir=pdb_dir, output_dir=output_dir)
    
    with multiprocessing.Pool(processes=N) as pool:
        for result in tqdm(pool.imap_unordered(process_func, pdb_list), total=len(pdb_list)):
            status, pdbid, data = result
            if status == 'success':
                log.writelines(data)
            else:
                log.write(f"{pdbid} {data}\n")
            log.flush()
    
    log.close()
    