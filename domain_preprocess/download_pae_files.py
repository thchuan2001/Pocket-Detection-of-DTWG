#!/usr/bin/env python3
"""
下载所有 AlphaFold 蛋白质的 PAE (Predicted Aligned Error) JSON 文件
"""

import os
import requests
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import time

def extract_protein_id_from_filename(filename):
    """从PDB文件名提取蛋白质ID"""
    # 文件名格式: AF-A0A385XJ53-F1-model_v6.pdb
    if filename.startswith('AF-') and '-F1-' in filename:
        parts = filename.split('-')
        if len(parts) >= 3:
            return parts[1]  # 返回 UniProt ID
    return None

def download_pae_file(args):
    """下载单个 PAE JSON 文件"""
    protein_id, output_dir, max_retries = args
    
    # 构建URL和输出文件路径
    url = f"https://alphafold.ebi.ac.uk/files/AF-{protein_id}-F1-predicted_aligned_error_v6.json"
    output_file = os.path.join(output_dir, f"AF-{protein_id}-F1-predicted_aligned_error_v6.json")
    
    # 如果文件已存在，跳过
    if os.path.exists(output_file):
        return (protein_id, 'skipped', None)
    
    # 尝试下载
    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            
            if response.status_code == 200:
                # 保存文件
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                return (protein_id, 'success', None)
            elif response.status_code == 404:
                return (protein_id, 'not_found', 'File not found on server')
            else:
                if attempt < max_retries - 1:
                    time.sleep(1)  # 重试前等待
                    continue
                return (protein_id, 'failed', f'HTTP {response.status_code}')
                
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(1)
                continue
            return (protein_id, 'error', str(e))
    
    return (protein_id, 'failed', 'Max retries reached')

def collect_protein_ids(base_dir):
    """收集所有蛋白质ID"""
    base_path = Path(base_dir)
    protein_ids = set()
    
    print("正在扫描PDB文件...")
    for subdir in base_path.iterdir():
        if not subdir.is_dir():
            continue
        
        print(f"  扫描 {subdir.name}...")
        for file_path in subdir.iterdir():
            if file_path.name.endswith('.pdb.gz') or file_path.name.endswith('.pdb'):
                protein_id = extract_protein_id_from_filename(file_path.name)
                if protein_id:
                    protein_ids.add((protein_id, subdir.name))
    
    return protein_ids

def main():
    base_dir = 'af2db_downloads'
    max_retries = 3
    
    if not Path(base_dir).exists():
        print(f"错误: 目录 {base_dir} 不存在")
        return
    
    print("=" * 80)
    print("AlphaFold PAE JSON 文件下载")
    print("=" * 80)
    
    # 收集所有蛋白质ID
    protein_data = collect_protein_ids(base_dir)
    
    print(f"\n找到 {len(protein_data)} 个唯一的蛋白质ID")
    
    # 按物种组织输出目录
    species_map = {}
    for protein_id, species in protein_data:
        if species not in species_map:
            species_map[species] = []
        species_map[species].append(protein_id)
    
    print(f"\n各物种分布:")
    for species, ids in sorted(species_map.items()):
        print(f"  {species}: {len(ids)} 个蛋白质")
    
    # 为每个物种创建PAE输出目录
    for species in species_map.keys():
        pae_dir = os.path.join(base_dir, species, 'pae_files')
        os.makedirs(pae_dir, exist_ok=True)
    
    # 准备下载参数
    download_args = []
    for protein_id, species in protein_data:
        output_dir = os.path.join(base_dir, species, 'pae_files')
        download_args.append((protein_id, output_dir, max_retries))
    
    # 确认操作
    print(f"\n将下载 {len(download_args)} 个 PAE JSON 文件")
    response = input("确认继续? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("操作已取消")
        return
    
    # 使用多进程下载
    n_processes = min(cpu_count(), 16)  # 限制并发数避免服务器压力
    print(f"\n使用 {n_processes} 个进程下载...")
    
    success_count = 0
    skipped_count = 0
    failed_count = 0
    not_found_count = 0
    error_count = 0
    
    with Pool(n_processes) as pool:
        with tqdm(total=len(download_args), desc="下载 PAE") as pbar:
            for protein_id, status, error in pool.imap_unordered(download_pae_file, download_args):
                if status == 'success':
                    success_count += 1
                elif status == 'skipped':
                    skipped_count += 1
                elif status == 'not_found':
                    not_found_count += 1
                elif status == 'failed':
                    failed_count += 1
                elif status == 'error':
                    error_count += 1
                
                pbar.update(1)
                pbar.set_postfix({
                    'success': success_count,
                    'skipped': skipped_count,
                    'not_found': not_found_count,
                    'failed': failed_count + error_count
                })
    
    # 输出统计
    print("\n" + "=" * 80)
    print("下载完成!")
    print("=" * 80)
    print(f"  成功下载: {success_count}")
    print(f"  已存在跳过: {skipped_count}")
    print(f"  服务器未找到: {not_found_count}")
    print(f"  下载失败: {failed_count}")
    print(f"  网络错误: {error_count}")
    print("=" * 80)

if __name__ == '__main__':
    main()
