#!/usr/bin/env python3
"""
对 Human 蛋白质运行 domain_split (无质量过滤)
"""

import os
import sys

# 导入 domain_split 模块
sys.path.insert(0, '/home/tanhaichuan/GenPack_for_galaxy')
from domain_split import main

if __name__ == '__main__':
    base_dir = '/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads'
    human_dir = os.path.join(base_dir, 'Human')
    
    # 生成 Human PDB 文件列表
    pdb_list_file = os.path.join(base_dir, 'human_pdb_list.txt')
    pdb_files = [f for f in os.listdir(human_dir) if f.endswith('.pdb')]
    
    print(f"找到 {len(pdb_files)} 个 Human PDB 文件")
    
    with open(pdb_list_file, 'w') as f:
        for pdb_file in sorted(pdb_files):
            f.write(f"{pdb_file}\n")
    
    print(f"PDB 文件列表已保存: {pdb_list_file}")
    
    # 配置路径
    protein_id_list = pdb_list_file
    pdb_dir = human_dir
    pae_dir = os.path.join(base_dir, 'Human_pae_files')  # PAE 文件需要单独下载
    output_dir = os.path.join(base_dir, 'Human_domains')
    N_CPU = 32
    total_files = len(pdb_files)
    log_path = os.path.join(base_dir, 'human_domain_split_errors.txt')
    
    # 检查 PAE 目录
    if not os.path.exists(pae_dir):
        print(f"\n警告: PAE 目录不存在: {pae_dir}")
        print("需要先下载 Human 的 PAE 文件")
        sys.exit(1)
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n配置信息:")
    print(f"  PDB 目录: {pdb_dir}")
    print(f"  PAE 目录: {pae_dir}")
    print(f"  输出目录: {output_dir}")
    print(f"  进程数: {N_CPU}")
    print(f"  总文件数: {total_files}")
    print(f"  错误日志: {log_path}")
    
    response = input("\n确认开始处理? (yes/no): ")
    if response.lower() not in ['yes', 'y']:
        print("操作已取消")
        sys.exit(0)
    
    print("\n开始处理...")
    main(protein_id_list, total_files, pdb_dir, pae_dir, output_dir, N_CPU, log_path)
    
    print("\n处理完成！")
