#!/usr/bin/env python3
"""
下载 Human 蛋白质的 PAE JSON 文件
"""

import os
import requests
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import time

def extract_protein_id(filename):
    """从PDB文件名提取蛋白质ID"""
    # AF-A0A385XJ53-F1-model_v6.pdb -> A0A385XJ53
    if filename.startswith('AF-') and '-F1-' in filename:
        parts = filename.split('-')
        if len(parts) >= 3:
            return parts[1]
    return None

def download_pae_file(args):
    """下载单个 PAE JSON 文件"""
    protein_id, output_dir, max_retries = args
    
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
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                return (protein_id, 'success', None)
            elif response.status_code == 404:
                return (protein_id, 'not_found', 'File not found on server')
            else:
                if attempt < max_retries - 1:
                    time.sleep(1)
                    continue
                return (protein_id, 'failed', f'HTTP {response.status_code}')
                
        except Exception as e:
            if attempt < max_retries - 1:
                time.sleep(1)
                continue
            return (protein_id, 'error', str(e))
    
    return (protein_id, 'failed', 'Max retries reached')

def main():
    base_dir = Path('/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads')
    human_pdb_dir = base_dir / 'Human'
    output_dir = base_dir / 'Human_pae_files'
    
    print("=" * 80)
    print("下载 Human PAE JSON 文件")
    print("=" * 80)
    
    # 创建输出目录
    output_dir.mkdir(exist_ok=True)
    
    # 收集所有蛋白质ID
    print("\n扫描 Human PDB 文件...")
    protein_ids = []
    for filename in os.listdir(human_pdb_dir):
        if filename.endswith('.pdb'):
            protein_id = extract_protein_id(filename)
            if protein_id:
                protein_ids.append(protein_id)
    
    print(f"找到 {len(protein_ids)} 个蛋白质ID")
    
    # 准备下载参数
    max_retries = 3
    download_args = [(pid, str(output_dir), max_retries) for pid in protein_ids]
    
    response = input(f"\n确认下载 {len(download_args)} 个 PAE 文件? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("操作已取消")
        return
    
    # 使用多进程下载
    n_processes = min(cpu_count(), 16)
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
