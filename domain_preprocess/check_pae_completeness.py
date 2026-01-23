#!/usr/bin/env python3
"""
检查 PDB 和 PAE 文件的配对情况
"""

import os
from pathlib import Path

def extract_protein_id(filename):
    """从文件名提取蛋白质ID"""
    # AF-A0A385XJ53-F1-model_v6.pdb -> A0A385XJ53
    if filename.startswith('AF-') and '-F1-' in filename:
        parts = filename.split('-')
        if len(parts) >= 3:
            return parts[1]
    return None

def main():
    base_dir = Path('/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads')
    pdb_dir = base_dir / 'pdb_files'
    pae_dir = base_dir / 'pae_files'
    
    print("正在扫描文件...")
    
    # 收集PDB文件的蛋白质ID
    pdb_ids = set()
    for filename in os.listdir(pdb_dir):
        if filename.endswith('.pdb'):
            protein_id = extract_protein_id(filename)
            if protein_id:
                pdb_ids.add(protein_id)
    
    # 收集PAE文件的蛋白质ID
    pae_ids = set()
    for filename in os.listdir(pae_dir):
        if filename.endswith('.json'):
            protein_id = extract_protein_id(filename)
            if protein_id:
                pae_ids.add(protein_id)
    
    # 分析差异
    missing_pae = pdb_ids - pae_ids
    extra_pae = pae_ids - pdb_ids
    matched = pdb_ids & pae_ids
    
    print("\n" + "=" * 80)
    print("文件配对检查结果")
    print("=" * 80)
    print(f"PDB 文件数量:     {len(pdb_ids):>10}")
    print(f"PAE 文件数量:     {len(pae_ids):>10}")
    print(f"配对成功:         {len(matched):>10}")
    print(f"缺少 PAE 文件:   {len(missing_pae):>10}")
    print(f"多余 PAE 文件:   {len(extra_pae):>10}")
    print("=" * 80)
    
    if missing_pae:
        print(f"\n缺少 PAE 文件的蛋白质 ID (前20个):")
        for i, protein_id in enumerate(sorted(missing_pae)[:20], 1):
            print(f"  {i:2}. {protein_id}")
        
        if len(missing_pae) > 20:
            print(f"  ... 还有 {len(missing_pae) - 20} 个")
        
        # 保存完整列表
        missing_file = base_dir / 'missing_pae_ids.txt'
        with open(missing_file, 'w') as f:
            for protein_id in sorted(missing_pae):
                f.write(f"{protein_id}\n")
        print(f"\n完整列表已保存到: {missing_file}")
    
    if extra_pae:
        print(f"\n多余 PAE 文件 (前20个):")
        for i, protein_id in enumerate(sorted(extra_pae)[:20], 1):
            print(f"  {i:2}. {protein_id}")
    
    print("\n" + "=" * 80)
    if len(missing_pae) == 0 and len(extra_pae) == 0:
        print("✓ 所有 PDB 文件都有对应的 PAE 文件")
    else:
        print("✗ 存在不匹配的文件")
    print("=" * 80)

if __name__ == '__main__':
    main()
