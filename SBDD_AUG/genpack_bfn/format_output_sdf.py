import os
import glob
import rdkit 
from rdkit import Chem
from tqdm import tqdm
import random
import shutil
import multiprocessing

def process_case(case_name, ligands, output_exps, max_ligands):
    output_dir = os.path.join(output_exps, case_name)
    os.makedirs(output_dir, exist_ok=True)
    if max_ligands:
        ligands = ligands[:max_ligands]
    for i, ligand in enumerate(ligands):
        cmd = 'cp ' + ligand + ' ' + os.path.join(output_dir, case_name + '_' + str(i) + '.sdf')
        os.system(cmd)

def format_output(
    input_exps= '/data/bfn_data/Flex_bfn/CrossDocked_backbone_APObind/default/test_outputs_v1/good_cases',
    output_exps= '/data/bfn_data/CrossDocked_sample_output/CrossDocked_backbone_APObind',
    max_ligands = None
):

    # if os.path.exists(output_exps):
    #     shutil.rmtree(output_exps)
    os.makedirs(output_exps, exist_ok=True)
    sdf_list = {}
    samples = glob.glob(input_exps + '/*')

    for sample in samples:

        # check one fragment only
        try:
            rdmol=Chem.MolFromMolFile(sample)
            assert rdmol is not None
            num_frags=Chem.GetMolFrags(rdmol,asMols=True)
            if len(num_frags)>1:
                raise ValueError('More than one fragment')
            
            case_name = "_".join(sample.split('/')[-1].split('_')[:-1]).replace(".sdf", "")
            if case_name not in sdf_list:
                sdf_list[case_name] = []
            sdf_list[case_name].append(sample)
        except Exception as e:
            print(f'Ignoring failure {sample}: {e}')
            continue

    with multiprocessing.Pool() as pool:
        pool.starmap(process_case, [(case_name, ligands, output_exps,max_ligands) for case_name, ligands in sdf_list.items()])

if __name__ == '__main__':
    format_output()