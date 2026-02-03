import os
import sys
from weakref import ref
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool
from Bio.PDB import PDBParser, Selection, PDBIO
from Bio.PDB.Polypeptide import is_aa
import glob
from copy import deepcopy
import subprocess
import pickle
import json

class TMaligner:
    def __init__(self):
        pass
    
    def run_TMalign(self, protein_file, ref_protein_file):

        rotation_matrix_file = os.path.join("/tmp",protein_file.split("/")[-1].split(".")[0]+"_"+ref_protein_file.split("/")[-1].split(".")[0]+"_rotation_matrix.txt")

        out_bytes = subprocess.check_output(['/data/template_matching_scripts/TMalign',protein_file,ref_protein_file,"-m",rotation_matrix_file])
        out_text = out_bytes.decode('utf-8').strip().split("\n")
        TMscore1=float(out_text[12].split(" ")[1])
        TMscore2=float(out_text[13].split(" ")[1])

        seq_protein = out_text[17]
        seq_ref_protein = out_text[19]

        return TMscore1,TMscore2,rotation_matrix_file ,seq_protein, seq_ref_protein


    def get_rotate_matrix(self,rotate_matrix_file):
        with open(rotate_matrix_file,"r") as f:
            data=f.readlines()
        u=[]
        t=[]
        for i in range(2,5):
            line=data[i].split(" ")
            line_float=[float(x) for x in line if x!=""]
            t.append(line_float[1])
            u.append(line_float[2:])
        u=np.array(u)
        t=np.array(t)
        # rm file
        os.remove(rotate_matrix_file)
        return u,t


    def rotate_protein(self, protein_file, rotation_matrix_file, output_file):
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", protein_file)
        model = structure[0]
        
        # read rotate_matrix
        rotation_matrix=self.get_rotate_matrix(rotation_matrix_file)
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coord=atom.get_coord()
                    coord=np.array(coord)
                    new_coord=np.dot(rotation_matrix[0],coord)+rotation_matrix[1]
                    atom.set_coord(new_coord)
        
        # write new pdb file
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_file)

    def align(self,protein_file,ref_protein_file):
        TMscore1,TMscore2,rotation_matrix_file,seq_protein,seq_ref_protein =self.run_TMalign(protein_file,ref_protein_file)
        rotation_matrix=self.get_rotate_matrix(rotation_matrix_file)
        return {
            "TMscore1":TMscore1,
            "TMscore2":TMscore2,
            "rotation_matrix":(rotation_matrix[0].tolist(), rotation_matrix[1].tolist()),
            "seq_protein":seq_protein,
            "seq_ref_protein":seq_ref_protein
        }
    

if __name__ == "__main__":

    # for manual multi node calculation
    thread_cnt=20
    thread_id=0

    # AF2_dir is the directory of AF2 domains. Domains are stored in pdb format.
    # PDBBind_list is the directory of PDBBind proteins. Proteins are stored in pdb format.

    AF2_dir="/data/domains"
    PDBBind_list="/data/DTWG_pdbbind_receptor_only"
    output_dir="/data/TMalign_results"

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    AF2_list=glob.glob(AF2_dir+"/*.pdb")
    PDBBind_list=glob.glob(PDBBind_list+"/*.pdb")
    AF2_list.sort()
    PDBBind_list.sort()

    # for manual multi node calculation
    new_AF2_list=[]
    for i,AF2 in enumerate(AF2_list):
        if i%thread_cnt==thread_id:
            new_AF2_list.append(AF2)
    AF2_list=new_AF2_list


    print("number of AF2 proteins:",len(AF2_list))
    print("number of PDBBind proteins:",len(PDBBind_list))

    def run(AF2_item):
        res=[]
        tmaligner=TMaligner()
        for pdbbind_item in tqdm(PDBBind_list):
            result=tmaligner.align(pdbbind_item,AF2_item)
            result["AF2"]=AF2_item.split("/")[-1].split(".")[0]
            result["PDBBind"]=pdbbind_item.split("/")[-1].split(".")[0]
            if result["TMscore2"]>0.5:
                res.append(result)
        # save 
        # with open(os.path.join(output_dir,AF2_item.split("/")[-1].split(".")[0]+".pkl"),"wb") as f:
        #     pickle.dump(res,f)

        # save as json
        with open(os.path.join(output_dir,AF2_item.split("/")[-1].split(".")[0]+".json"),"w") as f:
            json.dump(res,f)

    
    with Pool(500) as p:
        p.map(run, AF2_list)

    print("done")