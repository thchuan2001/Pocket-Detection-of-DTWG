#!/usr/bin/env python3
"""
处理 Human 文件夹：删除所有 .cif.gz 文件，解压所有 .pdb.gz 文件
"""

import os
import gzip
import shutil
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

def delete_cif_gz(file_path):
    """删除单个 .cif.gz 文件"""
    try:
        os.remove(file_path)
        return (str(file_path), 'deleted', None)
    except Exception as e:
        return (str(file_path), 'failed', str(e))

def decompress_pdb_gz(file_path):
    """解压单个 .pdb.gz 文件并删除原文件"""
    try:
        output_path = str(file_path).replace('.pdb.gz', '.pdb')
        
        with gzip.open(file_path, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # 删除原始 .gz 文件
        os.remove(file_path)
        
        return (str(file_path), 'decompressed', None)
    except Exception as e:
        return (str(file_path), 'failed', str(e))

def collect_files(human_dir):
    """收集所有需要处理的文件"""
    cif_gz_files = []
    pdb_gz_files = []
    
    print("正在扫描 Human 目录...")
    for root, dirs, files in os.walk(human_dir):
        for filename in files:
            file_path = Path(root) / filename
            if filename.endswith('.cif.gz'):
                cif_gz_files.append(file_path)
            elif filename.endswith('.pdb.gz'):
                pdb_gz_files.append(file_path)
    
    return cif_gz_files, pdb_gz_files

def main():
    human_dir = Path('/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads/Human')
    
    if not human_dir.exists():
        print(f"错误: 目录 {human_dir} 不存在")
        return
    
    print("=" * 80)
    print("处理 Human 文件夹")
    print("=" * 80)
    
    # 收集文件
    cif_gz_files, pdb_gz_files = collect_files(human_dir)
    
    print(f"\n找到 {len(cif_gz_files)} 个 .cif.gz 文件")
    print(f"找到 {len(pdb_gz_files)} 个 .pdb.gz 文件")
    
    if len(cif_gz_files) == 0 and len(pdb_gz_files) == 0:
        print("\n没有需要处理的文件")
        return
    
    print(f"\n将要执行的操作:")
    print(f"  1. 删除 {len(cif_gz_files)} 个 .cif.gz 文件")
    print(f"  2. 解压 {len(pdb_gz_files)} 个 .pdb.gz 文件")
    print(f"  3. 删除原始 .pdb.gz 文件")
    
    response = input("\n确认继续? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("操作已取消")
        return
    
    n_processes = min(cpu_count(), 32)
    print(f"\n使用 {n_processes} 个进程...")
    
    # 删除 .cif.gz 文件
    if cif_gz_files:
        print(f"\n正在删除 {len(cif_gz_files)} 个 .cif.gz 文件...")
        with Pool(n_processes) as pool:
            results = list(tqdm(
                pool.imap(delete_cif_gz, cif_gz_files),
                total=len(cif_gz_files),
                desc="删除 CIF.gz"
            ))
        
        deleted = sum(1 for _, status, _ in results if status == 'deleted')
        failed = sum(1 for _, status, _ in results if status == 'failed')
        print(f"  成功删除: {deleted}")
        if failed > 0:
            print(f"  失败: {failed}")
    
    # 解压 .pdb.gz 文件
    if pdb_gz_files:
        print(f"\n正在解压 {len(pdb_gz_files)} 个 .pdb.gz 文件...")
        with Pool(n_processes) as pool:
            results = list(tqdm(
                pool.imap(decompress_pdb_gz, pdb_gz_files),
                total=len(pdb_gz_files),
                desc="解压 PDB.gz"
            ))
        
        decompressed = sum(1 for _, status, _ in results if status == 'decompressed')
        failed = sum(1 for _, status, _ in results if status == 'failed')
        print(f"  成功解压: {decompressed}")
        if failed > 0:
            print(f"  失败: {failed}")
    
    print("\n" + "=" * 80)
    print("处理完成！")
    print("=" * 80)

if __name__ == '__main__':
    main()
