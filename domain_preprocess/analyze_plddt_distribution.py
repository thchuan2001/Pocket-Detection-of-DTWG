#!/usr/bin/env python3
"""
统计不同物种 AlphaFold 蛋白质的 pLDDT 分布
"""

import numpy as np
from Bio.PDB import PDBParser
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import matplotlib.pyplot as plt
import pandas as pd

def calculate_plddt_stats(pdb_path):
    """计算单个PDB文件的pLDDT统计信息"""
    try:
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', pdb_path)
        chain = structure[0]['A']
        
        plddt = []
        for residue in chain:
            if 'CA' in residue:
                plddt.append(residue['CA'].get_bfactor())
        
        plddt = np.array(plddt)
        
        if len(plddt) == 0:
            return None
        
        mean_plddt = plddt.mean()
        ratio_above_70 = sum(plddt >= 70) / len(plddt)
        
        return {
            'filename': os.path.basename(pdb_path),
            'mean_plddt': mean_plddt,
            'ratio_above_70': ratio_above_70,
            'num_residues': len(plddt)
        }
    except Exception as e:
        return None

def collect_pdb_files_from_directory(directory):
    """收集目录下所有PDB文件"""
    pdb_files = []
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith('.pdb'):
                pdb_files.append(os.path.join(root, filename))
    return pdb_files

def collect_pdb_files_from_id_list(id_list_file, pdb_dir):
    """根据ID列表文件收集PDB文件路径"""
    pdb_files = []
    with open(id_list_file, 'r') as f:
        for line in f:
            protein_id = line.strip()
            if protein_id:
                # 构造PDB文件名: AF-{ID}-F1-model_v6.pdb
                pdb_filename = f"AF-{protein_id}-F1-model_v6.pdb"
                pdb_path = os.path.join(pdb_dir, pdb_filename)
                if os.path.exists(pdb_path):
                    pdb_files.append(pdb_path)
    return pdb_files

def analyze_species(species_name, pdb_source, n_processes=32, use_id_list=False, pdb_dir=None):
    """分析单个物种的pLDDT分布
    
    Args:
        species_name: 物种名称
        pdb_source: PDB文件来源（目录路径或ID列表文件路径）
        n_processes: 进程数
        use_id_list: 是否使用ID列表文件
        pdb_dir: 当use_id_list=True时，PDB文件所在的目录
    """
    print(f"\n{'='*80}")
    print(f"分析 {species_name}")
    print('='*80)
    
    # 收集PDB文件
    if use_id_list:
        pdb_files = collect_pdb_files_from_id_list(pdb_source, pdb_dir)
    else:
        pdb_files = collect_pdb_files_from_directory(pdb_source)
    
    print(f"找到 {len(pdb_files)} 个 PDB 文件")
    
    if len(pdb_files) == 0:
        return None
    
    # 并行处理
    print(f"使用 {n_processes} 个进程分析...")
    with Pool(n_processes) as pool:
        results = list(tqdm(
            pool.imap(calculate_plddt_stats, pdb_files),
            total=len(pdb_files),
            desc=f"分析 {species_name}"
        ))
    
    # 过滤掉失败的结果
    results = [r for r in results if r is not None]
    
    if len(results) == 0:
        print("没有成功分析的文件")
        return None
    
    # 转换为DataFrame
    df = pd.DataFrame(results)
    
    # 统计信息
    print(f"\n统计结果:")
    print(f"  成功分析: {len(df)} 个蛋白")
    print(f"  平均 pLDDT 均值: {df['mean_plddt'].mean():.2f} ± {df['mean_plddt'].std():.2f}")
    print(f"  平均 pLDDT 中位数: {df['mean_plddt'].median():.2f}")
    print(f"  平均 ratio≥70: {df['ratio_above_70'].mean():.4f} ± {df['ratio_above_70'].std():.4f}")
    print(f"  平均 ratio≥70 中位数: {df['ratio_above_70'].median():.4f}")
    
    # 质量过滤统计
    high_quality = df[(df['mean_plddt'] >= 69.9) & (df['ratio_above_70'] >= 0.899)]
    print(f"\n质量过滤 (mean_pLDDT≥69.9 且 ratio≥70≥0.899):")
    print(f"  通过过滤: {len(high_quality)} ({len(high_quality)/len(df)*100:.2f}%)")
    print(f"  被过滤: {len(df) - len(high_quality)} ({(len(df)-len(high_quality))/len(df)*100:.2f}%)")
    
    return df

def plot_distributions(species_data, output_dir):
    """绘制分布图"""
    print(f"\n{'='*80}")
    print("生成分布图...")
    print('='*80)
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 设置图表样式
    plt.style.use('default')
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    # 1. 平均 pLDDT 分布
    fig, axes = plt.subplots(len(species_data), 2, figsize=(15, 5*len(species_data)))
    if len(species_data) == 1:
        axes = axes.reshape(1, -1)
    
    for idx, (species_name, df) in enumerate(species_data.items()):
        color = colors[idx % len(colors)]
        
        # Mean pLDDT distribution
        axes[idx, 0].hist(df['mean_plddt'], bins=50, color=color, alpha=0.7, edgecolor='black')
        axes[idx, 0].axvline(69.9, color='red', linestyle='--', linewidth=2, label='Threshold: 69.9')
        axes[idx, 0].set_xlabel('Mean pLDDT', fontsize=12)
        axes[idx, 0].set_ylabel('Number of Proteins', fontsize=12)
        axes[idx, 0].set_title(f'{species_name} - Mean pLDDT Distribution', fontsize=14, fontweight='bold')
        axes[idx, 0].legend()
        axes[idx, 0].grid(alpha=0.3)
        
        # ratio≥70 distribution
        axes[idx, 1].hist(df['ratio_above_70'], bins=50, color=color, alpha=0.7, edgecolor='black')
        axes[idx, 1].axvline(0.899, color='red', linestyle='--', linewidth=2, label='Threshold: 0.899')
        axes[idx, 1].set_xlabel('Ratio of Residues with pLDDT>=70', fontsize=12)
        axes[idx, 1].set_ylabel('Number of Proteins', fontsize=12)
        axes[idx, 1].set_title(f'{species_name} - pLDDT>=70 Ratio Distribution', fontsize=14, fontweight='bold')
        axes[idx, 1].legend()
        axes[idx, 1].grid(alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'plddt_distributions.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"保存图表: {output_file}")
    plt.close()
    
    # 2. 对比图（如果有多个物种）
    if len(species_data) > 1:
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))
        
        for idx, (species_name, df) in enumerate(species_data.items()):
            color = colors[idx % len(colors)]
            axes[0].hist(df['mean_plddt'], bins=50, alpha=0.5, label=species_name, color=color)
            axes[1].hist(df['ratio_above_70'], bins=50, alpha=0.5, label=species_name, color=color)
        
        axes[0].axvline(69.9, color='red', linestyle='--', linewidth=2, label='Threshold: 69.9')
        axes[0].set_xlabel('Mean pLDDT', fontsize=12)
        axes[0].set_ylabel('Number of Proteins', fontsize=12)
        axes[0].set_title('All Species - Mean pLDDT Distribution Comparison', fontsize=14, fontweight='bold')
        axes[0].legend()
        axes[0].grid(alpha=0.3)
        
        axes[1].axvline(0.899, color='red', linestyle='--', linewidth=2, label='Threshold: 0.899')
        axes[1].set_xlabel('Ratio of Residues with pLDDT>=70', fontsize=12)
        axes[1].set_ylabel('Number of Proteins', fontsize=12)
        axes[1].set_title('All Species - pLDDT>=70 Ratio Distribution Comparison', fontsize=14, fontweight='bold')
        axes[1].legend()
        axes[1].grid(alpha=0.3)
        
        plt.tight_layout()
        output_file = os.path.join(output_dir, 'plddt_distributions_comparison.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"保存对比图: {output_file}")
        plt.close()

def main():
    base_dir = Path('/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads')
    output_dir = base_dir / 'plddt_analysis'
    species_lists_dir = base_dir / 'species_lists'
    pdb_files_dir = base_dir / 'pdb_files'
    
    # 定义要分析的物种
    # Human单独存储在自己的目录中
    # 其他5个物种的PDB文件都在pdb_files目录，通过species_lists中的ID列表区分
    species_config = {
        'Human': {
            'type': 'directory',
            'source': base_dir / 'Human'
        }
    }
    
    # 添加其他5个物种（通过ID列表）
    if species_lists_dir.exists():
        for id_file in species_lists_dir.glob('*_ids.txt'):
            species_name = id_file.stem.replace('_ids', '')
            species_config[species_name] = {
                'type': 'id_list',
                'source': id_file,
                'pdb_dir': pdb_files_dir
            }
    
    # 检查并准备分析
    species_to_analyze = []
    for species_name, config in species_config.items():
        if config['type'] == 'directory':
            if config['source'].exists():
                species_to_analyze.append((species_name, config))
                print(f"✓ 找到 {species_name}: {config['source']}")
            else:
                print(f"✗ 未找到 {species_name}: {config['source']}")
        else:  # id_list
            if config['source'].exists() and config['pdb_dir'].exists():
                species_to_analyze.append((species_name, config))
                print(f"✓ 找到 {species_name}: {config['source']}")
            else:
                print(f"✗ 未找到 {species_name}: {config['source']}")
    
    if not species_to_analyze:
        print("\n错误: 没有找到可分析的物种数据")
        return
    
    # 分析每个物种
    n_processes = min(cpu_count(), 32)
    species_data = {}
    
    for species_name, config in species_to_analyze:
        if config['type'] == 'directory':
            df = analyze_species(species_name, config['source'], n_processes, use_id_list=False)
        else:  # id_list
            df = analyze_species(species_name, config['source'], n_processes, use_id_list=True, pdb_dir=config['pdb_dir'])
        
        if df is not None:
            species_data[species_name] = df
            # 保存详细数据
            os.makedirs(output_dir, exist_ok=True)
            csv_file = output_dir / f'{species_name}_plddt_stats.csv'
            df.to_csv(csv_file, index=False)
            print(f"保存数据: {csv_file}")
    
    # 绘制分布图
    if species_data:
        plot_distributions(species_data, output_dir)
    
    print(f"\n{'='*80}")
    print("分析完成！")
    print(f"结果保存在: {output_dir}")
    print('='*80)

if __name__ == '__main__':
    main()
