import os
import sys
import glob
import pickle
import rdkit
from rdkit import Chem
import numpy as np
import json
import Bio
from Bio.PDB import PDBParser, Selection, PDBIO
from tqdm import tqdm
from multiprocessing import Pool
import random

AA_NAME_SYM = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H',
    'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q',
    'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}

def remove_gaps(seq1,seq2):
    assert len(seq1)==len(seq2)
    new_seq1=""
    new_seq2=""
    for i in range(len(seq1)):
        if seq1[i]!="-" :
            new_seq1+=seq1[i]
            new_seq2+=seq2[i]
    return new_seq1,new_seq2

def get_pocket_match_rate(item,AF2_dir,PDBBind_dir,ligand_dir,pocket_position_dir):

    pdbbind_seq=item['seq_protein']
    AF2_seq=item['seq_ref_protein']

    pdbbbind_seq,AF2_seq=remove_gaps(pdbbind_seq,AF2_seq)
    position_file=os.path.join(pocket_position_dir,item["PDBBind"]+".txt")
    if not os.path.exists(position_file):
        raise Exception("position file not found")
    with open(position_file,"r") as f:
        pdbbbind_seq2=f.readline().strip()
        pocket_pos=f.readline().strip()
    
    assert pdbbbind_seq==pdbbbind_seq2

    match_cnt=0
    total_cnt=0
    for i in range(len(pdbbbind_seq)):
        if pocket_pos[i]!="-":
            total_cnt+=1
            if AF2_seq[i]!="-":
                match_cnt+=1

    return match_cnt/total_cnt,match_cnt,total_cnt


def rotate_ligand(rd_mol, rotation_matrix):
    conf = rd_mol.GetConformer()
    ligand_coords = conf.GetPositions()

    for i in range(len(ligand_coords)):
        ligand_coords[i] = np.dot(rotation_matrix[0], ligand_coords[i]) + rotation_matrix[1]

    for i in range(len(ligand_coords)):
        conf.SetAtomPosition(i, ligand_coords[i])

    return rd_mol


def get_pocket_ids(protein_file,ligand):
    ret=set()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_file)
    model = structure[0]
    conf = ligand.GetConformer()
    ligand_coords = conf.GetPositions()
    for chain in model:
        for residue in chain:
            f=0
            for atom in residue:
                protein_atom_coords = np.array(atom.coord)
                for ligand_coord in ligand_coords:
                    if np.linalg.norm(protein_atom_coords - ligand_coord) < 6.0:
                        f=1
                        break
                if f==1:
                    break
            if f:
                ret.add(residue.get_id()[1])
    return ret

def calc_iou(pocket_ids1,pocket_ids2):
    return len(pocket_ids1.intersection(pocket_ids2))/len(pocket_ids1.union(pocket_ids2))

def get_pos_to_id(protein_file):
    ret={}
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", protein_file)
    model = structure[0]
    cnt=0
    for chain in model:
        for residue in chain:
            ret[cnt]=residue.get_id()[1]
            cnt+=1
    return ret

def get_AF2_matched_ids(AF2_seq,pdbbind_seq,pos_to_id):
    AF2_seq,pdbbind_seq=remove_gaps(AF2_seq,pdbbind_seq)
    assert len(AF2_seq)==len(pos_to_id)
    ret=set()
    for i in range(len(pdbbind_seq)):
        if pdbbind_seq[i]!="-":
            ret.add(pos_to_id[i])
    return ret

def process_AF2_item(AF2_item, output_dir, new_TMalign_result_dir, ligand_dir, AF2_dir, PDBBind_dir, pocket_position_dir,single_chain_dir):
    new_res = []
    if AF2_item.endswith(".pkl"):
        with open(AF2_item, "rb") as f:
            res_list = pickle.load(f)
    elif AF2_item.endswith(".json"):
        with open(AF2_item, "r") as f:
            res_list = json.load(f)
    sample_size = min(1000, len(res_list))
    res_list = random.sample(res_list, sample_size)
    single_chain_ids=pickle.load(open(single_chain_dir,"rb"))
    pos_to_id=get_pos_to_id(os.path.join(AF2_dir, AF2_item.split("/")[-1].replace(".pkl",".pdb").replace(".json",".pdb")))
    pocket_cluster = []
    for item in tqdm(res_list):
        if "rotation_matrix" in item and isinstance(item["rotation_matrix"], list):
            item["rotation_matrix"] = (np.array(item["rotation_matrix"][0]), np.array(item["rotation_matrix"][1]))
            
        if item["PDBBind"] not in single_chain_ids:
            continue
        try:
            if item['TMscore2'] < 0.6:
                continue
            pocket_match_rate, match_cnt, pdbbind_cnt = get_pocket_match_rate(item, AF2_dir, PDBBind_dir, ligand_dir, pocket_position_dir)
            if pocket_match_rate == -1:
                continue
            if pocket_match_rate < 0.7:
                continue

            mol = Chem.MolFromMol2File(os.path.join(ligand_dir, item["PDBBind"] + ".mol2"))
            if mol is None:
                continue
            if mol.GetNumAtoms() > 800:
                continue
            mol = rotate_ligand(mol, item["rotation_matrix"])
            pocket_ids = get_pocket_ids(
                protein_file=os.path.join(AF2_dir, item["AF2"] + ".pdb"),
                ligand=mol
            )
            AF2_cnt = len(pocket_ids)
            AF2_matched_ids=get_AF2_matched_ids(item['seq_ref_protein'],item['seq_protein'],pos_to_id)
            match_cnt=len(pocket_ids.intersection(AF2_matched_ids))

            pocket_match_iou = match_cnt / (AF2_cnt + pdbbind_cnt - match_cnt)
            if pocket_match_iou < 0.6:
                continue

            new_clutser_flag = True
            for i in range(len(pocket_cluster)):
                iou = calc_iou(pocket_cluster[i]['pocket_ids'], pocket_ids)
                if iou > 0.4:
                    new_clutser_flag = False
                    if pocket_cluster[i]['cnt'] < 5:
                        mol.SetProp("_Name", "_".join([item["AF2"], item["PDBBind"], str(i), str(pocket_cluster[i]['cnt'])]))
                        w = Chem.SDWriter(os.path.join(output_dir, "_".join([item["AF2"], item["PDBBind"], str(i), str(pocket_cluster[i]['cnt'])]) + ".sdf"))
                        w.write(mol)
                        w.close()
                        pocket_cluster[i]['cnt'] += 1
                        print("found hit")
                        item["pocket_match_rate"] = pocket_match_rate
                        new_res.append(item)
                    else:
                        print("cluster full")
                    break
            if new_clutser_flag:
                mol.SetProp("_Name", "_".join([item["AF2"], item["PDBBind"], str(len(pocket_cluster)), "0"]))
                w = Chem.SDWriter(os.path.join(output_dir, "_".join([item["AF2"], item["PDBBind"], str(len(pocket_cluster)), "0"]) + ".sdf"))
                w.write(mol)
                w.close()
                pocket_cluster.append({
                    "pocket_ids": pocket_ids,
                    "cnt": 1
                })  
                print("found hit")
                item["pocket_match_rate"] = pocket_match_rate
                new_res.append(item)
                
        except KeyboardInterrupt:
            raise KeyboardInterrupt
        except Exception as e:
            print(e)
            pass

    new_result_file = os.path.join(new_TMalign_result_dir, os.path.basename(AF2_item))
    with open(new_result_file, "wb") as f:
        pickle.dump(new_res, f)
    return None

if __name__ == "__main__":
    TMalign_result_dir = "/data/TMalign_results"
    output_dir = "/data/pocket_match_results"
    new_TMalign_result_dir = "/data/pocket_match_filtered_results"
    ligand_dir = "/data/pdbbind_ligand_only"
    AF2_dir = "/data/domains"
    PDBBind_dir = "/data/DTWG_pdbbind_receptor_only"
    pocket_position_dir = "/data/pdbbind_pocket6A_position"
    single_chain_dir = "/data/single_chain_pocket10A.pkl"

    TMalign_results = glob.glob(TMalign_result_dir + "/*.json")

    # #######
    # # debug
    # TMalign_results = ["/data/DTWG_AF2_TM05/AF-P19835-F1-model_v4_0.pkl"]
    # #######


    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(new_TMalign_result_dir):
        os.makedirs(new_TMalign_result_dir)

    pool = Pool()

    pool.starmap(process_AF2_item, [(AF2_item, output_dir, new_TMalign_result_dir, ligand_dir, AF2_dir, PDBBind_dir, pocket_position_dir,single_chain_dir) for AF2_item in TMalign_results])

    ####### 
    # debug
    # for AF2_item in TMalign_results:
    #     process_AF2_item(AF2_item, output_dir, new_TMalign_result_dir, ligand_dir, AF2_dir, PDBBind_dir, pocket_position_dir,single_chain_dir)
    #######

    pool.close()
    pool.join()

    print("done")