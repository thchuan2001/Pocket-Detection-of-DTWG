#!/usr/bin/env python3
"""
重新组织 AlphaFold 数据库文件结构
- 记录每个物种的蛋白质ID列表
- 将所有PDB文件集中到 pdb_files/ 目录
- 将所有PAE文件集中到 pae_files/ 目录
"""

import os
import shutil
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm

def extract_protein_id(filename):
    """从文件名提取蛋白质ID"""
    # AF-A0A385XJ53-F1-model_v6.pdb.gz 或 AF-A0A385XJ53-F1-predicted_aligned_error_v6.json
    if filename.startswith('AF-') and '-F1-' in filename:
        parts = filename.split('-')
        if len(parts) >= 3:
            return parts[1]  # 返回 UniProt ID
    return None

def scan_current_structure(base_dir):
    """扫描当前目录结构，收集所有信息"""
    base_path = Path(base_dir)
    
    species_data = defaultdict(lambda: {
        'pdb_files': [],
        'pae_files': [],
        'protein_ids': set()
    })
    
    print("正在扫描现有文件结构...")
    
    for species_dir in base_path.iterdir():
        if not species_dir.is_dir():
            continue
        
        species_name = species_dir.name
        print(f"  扫描 {species_name}...")
        
        # 扫描主目录中的PDB文件
        for file_path in species_dir.iterdir():
            if file_path.is_file():
                filename = file_path.name
                protein_id = extract_protein_id(filename)
                
                if filename.endswith('.pdb') or filename.endswith('.pdb.gz'):
                    species_data[species_name]['pdb_files'].append(file_path)
                    if protein_id:
                        species_data[species_name]['protein_ids'].add(protein_id)
        
        # 扫描pae_files子目录
        pae_dir = species_dir / 'pae_files'
        if pae_dir.exists() and pae_dir.is_dir():
            for file_path in pae_dir.iterdir():
                if file_path.is_file() and file_path.name.endswith('.json'):
                    species_data[species_name]['pae_files'].append(file_path)
                    protein_id = extract_protein_id(file_path.name)
                    if protein_id:
                        species_data[species_name]['protein_ids'].add(protein_id)
    
    return species_data

def reorganize_files(base_dir, species_data):
    """重新组织文件结构"""
    base_path = Path(base_dir)
    
    # 创建新的目录结构
    pdb_dir = base_path / 'all_pdb_files'
    pae_dir = base_path / 'all_pae_files'
    species_lists_dir = base_path / 'species_lists'
    
    print("\n创建新目录结构...")
    pdb_dir.mkdir(exist_ok=True)
    pae_dir.mkdir(exist_ok=True)
    species_lists_dir.mkdir(exist_ok=True)
    
    # 写入每个物种的ID列表
    print("\n生成物种ID列表...")
    for species_name, data in sorted(species_data.items()):
        list_file = species_lists_dir / f"{species_name}_ids.txt"
        with open(list_file, 'w') as f:
            for protein_id in sorted(data['protein_ids']):
                f.write(f"{protein_id}\n")
        print(f"  {species_name}: {len(data['protein_ids'])} 个蛋白质ID")
    
    # 统计总文件数
    total_pdb = sum(len(data['pdb_files']) for data in species_data.values())
    total_pae = sum(len(data['pae_files']) for data in species_data.values())
    
    print(f"\n准备移动文件:")
    print(f"  PDB文件: {total_pdb} 个")
    print(f"  PAE文件: {total_pae} 个")
    
    # 移动PDB文件
    if total_pdb > 0:
        print(f"\n移动 {total_pdb} 个 PDB 文件...")
        with tqdm(total=total_pdb, desc="移动PDB") as pbar:
            for species_name, data in species_data.items():
                for pdb_file in data['pdb_files']:
                    dest = pdb_dir / pdb_file.name
                    if not dest.exists():
                        shutil.move(str(pdb_file), str(dest))
                    pbar.update(1)
    
    # 移动PAE文件
    if total_pae > 0:
        print(f"\n移动 {total_pae} 个 PAE 文件...")
        with tqdm(total=total_pae, desc="移动PAE") as pbar:
            for species_name, data in species_data.items():
                for pae_file in data['pae_files']:
                    dest = pae_dir / pae_file.name
                    if not dest.exists():
                        shutil.move(str(pae_file), str(dest))
                    pbar.update(1)
    
    # 清理空目录
    print("\n清理空目录...")
    for species_name in species_data.keys():
        species_dir = base_path / species_name
        pae_subdir = species_dir / 'pae_files'
        
        # 删除pae_files子目录（如果为空）
        if pae_subdir.exists():
            try:
                pae_subdir.rmdir()
                print(f"  删除空目录: {pae_subdir}")
            except OSError:
                print(f"  目录不为空，保留: {pae_subdir}")
        
        # 删除物种主目录（如果为空）
        try:
            species_dir.rmdir()
            print(f"  删除空目录: {species_dir}")
        except OSError:
            print(f"  目录不为空，保留: {species_dir}")
    
    return pdb_dir, pae_dir, species_lists_dir

def main():
    base_dir = 'af2db_downloads'
    
    if not Path(base_dir).exists():
        print(f"错误: 目录 {base_dir} 不存在")
        return
    
    print("=" * 80)
    print("AlphaFold 数据库文件结构重组")
    print("=" * 80)
    
    # 扫描当前结构
    species_data = scan_current_structure(base_dir)
    
    if not species_data:
        print("未找到任何数据文件")
        return
    
    # 显示统计信息
    print("\n当前文件统计:")
    for species_name, data in sorted(species_data.items()):
        print(f"  {species_name}:")
        print(f"    PDB文件: {len(data['pdb_files'])}")
        print(f"    PAE文件: {len(data['pae_files'])}")
        print(f"    唯一ID: {len(data['protein_ids'])}")
    
    # 确认操作
    print("\n" + "=" * 80)
    print("将执行以下操作:")
    print("  1. 在 species_lists/ 目录生成各物种的ID列表文件")
    print("  2. 将所有PDB文件移动到 all_pdb_files/ 目录")
    print("  3. 将所有PAE文件移动到 all_pae_files/ 目录")
    print("  4. 删除空的原始目录")
    print("=" * 80)
    
    response = input("\n确认继续? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("操作已取消")
        return
    
    # 重新组织文件
    pdb_dir, pae_dir, lists_dir = reorganize_files(base_dir, species_data)
    
    print("\n" + "=" * 80)
    print("重组完成!")
    print("=" * 80)
    print(f"PDB文件目录: {pdb_dir}")
    print(f"PAE文件目录: {pae_dir}")
    print(f"物种列表目录: {lists_dir}")
    print("=" * 80)

if __name__ == '__main__':
    main()
