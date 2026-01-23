# create a pkl file for bfn for test data

import os
import glob

# 使用10k测试口袋目录
fpocket_path="/home/tanhaichuan/GenPack_for_galaxy/test_10k_backbone_pocket"

# 输出到测试目录
pkl_output_path="/home/tanhaichuan/GenPack_for_galaxy/test_backbone_pocket.pkl"
ref_ch4_file="/project/ch4.sdf"
pockets=glob.glob(os.path.join(fpocket_path,"*.pdb"))

# 不使用关键词过滤，使用所有测试文件
print(f"Found {len(pockets)} test pockets")



result=[]

for pocket in pockets:
    pocket_name=os.path.basename(pocket).replace(".pdb","")
    ref_ligand_path=ref_ch4_file
    pocket_path=pocket

    result.append({
        "pocket_name":pocket_name,
        "pocket_path":pocket_path,
        "ref_ligand_path":ref_ligand_path
    })

import json
# print(json.dumps(result,indent=4))
import pickle
with open(pkl_output_path,"wb") as f:
    pickle.dump(result,f)
print("Finish")
