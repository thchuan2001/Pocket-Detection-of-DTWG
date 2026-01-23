from curses import raw
from core.evaluation.utils.file_utils import assemble_complex, disassemble_complex
from core.evaluation.utils.relax import relax
import concurrent.futures
from unittest import result
import torch
import rdkit
from rdkit import Chem
from tqdm import tqdm
import sys 
import os
import argparse
import glob
import subprocess
import pickle

"""新添加"""
import shutil


from eva_dataset import (
    TargetdiffCrossDockedCrossDockedUnrelaxed,
    TargetdiffCrossDockedBackboneCrossDockedUnrelaxed, 
    TargetdiffCrossDockedBackboneAPObindUnrelaxed, 
    TargetdiffCrossDockedAPObindUnrelaxed,
    BfnCrossDockedCrossDockedUnrelaxed,
    BfnCrossDockedBackboneCrossDockedUnrelaxed,
)
import os

class RelaxAndAssembleRunner:
    def __init__(self, dataset):
        self.dataset = dataset
        self.num_process = 100
        self.overwrite = True
    

    def process_single_item(self, item):
        for raw_sdf_path in tqdm(item["ligand_paths"]):
            try:
                print("processing:",raw_sdf_path)
                # make mol_path for each ligand
                mol_path = raw_sdf_path.replace(".sdf", "")
                if self.overwrite and os.path.exists(mol_path):
                    import shutil
                    shutil.rmtree(mol_path)
                os.makedirs(mol_path, exist_ok=True)

                # cp raw_sdf to mol_path
                os.system(f"cp {raw_sdf_path} {mol_path}")
                raw_sdf_path = os.path.join(mol_path, os.path.basename(raw_sdf_path))

                # cp protein to mol_path
                raw_protein_path = item["protein_path"]
                os.system(f"cp {raw_protein_path} {mol_path}")
                raw_protein_path = os.path.join(mol_path, os.path.basename(raw_protein_path))

                # assemble complex
                raw_complex_path = os.path.join(mol_path, "complex.pdb")
                if self.overwrite or not os.path.exists(raw_complex_path):
                    assemble_complex(raw_protein_path, raw_sdf_path, raw_complex_path)

                # relax
                refined_complex_path = raw_complex_path.replace(".pdb", "_refined.pdb")
                if self.overwrite or not os.path.exists(refined_complex_path):
                    relax(raw_complex_path)

                # disassemble complex
                refined_protein_path = raw_protein_path.replace(".pdb", "_refined.pdb")
                refined_sdf_path = raw_sdf_path.replace(".sdf", "_refined.sdf")

                if self.overwrite or not os.path.exists(refined_sdf_path):
                    disassemble_complex(refined_complex_path, refined_protein_path, refined_sdf_path)

            except KeyboardInterrupt:
                raise KeyboardInterrupt
            except Exception as e:
                print(f"###### Error: {raw_sdf_path}")
                print(e)
                pass
        
        print(f"Finish {item['name']}")

    def process_items(self, items):
        for item in items:
            self.process_single_item(item)
    
    def run(self):
        tasks = self.dataset.items

        # sort tasks by name
        if tasks[0]["name"].split("_")[0].isdigit():
            tasks = sorted(tasks, key=lambda x: int(x["name"].split("_")[0]))
        print(f"Total number of cases: {len(tasks)}")
        task_splits = [tasks[i::self.num_process] for i in range(self.num_process)]
        print("tasks_splits",task_splits)

        with concurrent.futures.ProcessPoolExecutor(max_workers=self.num_process) as executor:
            futures = []
            for i in range(self.num_process):
                future = executor.submit(self.process_items, task_splits[i])
                futures.append(future) 
            concurrent.futures.wait(futures)

        # for i in range(self.num_process):
        #     self.process_items(task_splits[i])
        print("Finish all cases")


class DiffDockandAssembleRunner:
    def __init__(self):
        self.gpu_ids=[1,2,3,4,5,6,7]
        self.overwrite = False
    
    def get_raw_protein_path(self,target_name):
        pkl=glob.glob(os.path.join(self.wd,"*.pkl"))
        if len(pkl)==0:
            print("Error: pkl not found")
            return None
        pkl=pkl[0]
        self.name2protein_path={}
        with open(pkl, "rb") as f:
            init_items = pickle.load(f)
            for item in init_items:
                self.name2protein_path[item["pocket_name"]]=item["protein_path"]
        return self.name2protein_path[target_name]


    def process_single_item(self, task, gpu_id):
        target_wd=task
        ligand_wds=[f for f in glob.glob(target_wd+"/*") if f.endswith(".sdf")]

        for raw_sdf_path in tqdm(ligand_wds):
            try:
                print("processing:",raw_sdf_path)
                # make mol_path for each ligand
                mol_path = raw_sdf_path.replace(".sdf", "")
                if self.overwrite and os.path.exists(mol_path):
                    import shutil
                    shutil.rmtree(mol_path)
                os.makedirs(mol_path, exist_ok=True)

                # cp raw_sdf to mol_path
                os.system(f"cp {raw_sdf_path} {mol_path}")
                raw_sdf_path = os.path.join(mol_path, os.path.basename(raw_sdf_path))

                # cp protein to mol_path
                source_protein_path=self.get_raw_protein_path(os.path.basename(target_wd))
                raw_protein_path = os.path.join(mol_path, os.path.basename(source_protein_path).replace("_refined",""))
                os.system(f"cp {source_protein_path} {raw_protein_path}")

                # skip if refined protein and sdf exist
                refined_protein_path = raw_protein_path.replace(".pdb", "_refined.pdb")
                refined_sdf_path = raw_sdf_path.replace(".sdf", "_refined.sdf")
                if os.path.exists(refined_protein_path) and os.path.exists(refined_sdf_path):
                    continue

                # run diffdock-pocket
                output_path=os.path.join(mol_path,"diffdock_output")
                python_path="/home/tanhaichuan/.conda/envs/diffdock_new_new/bin/python"
                cmd=f"cd /project/DiffDock-Pocket && CUDA_VISIBLE_DEVICES={gpu_id} {python_path} inference.py  --batch_size 10 --samples_per_complex 40 --keep_local_structures --model_dir /data/diffdock-pocket/model_ckpt/score_model --filtering_model_dir /data/diffdock-pocket/model_ckpt/confidence_model --protein_path {raw_protein_path} --ligand {raw_sdf_path} --out_dir {output_path} --relax"
                
                os.system(cmd)
                
                output_result=glob.glob(os.path.join(output_path,"*"))[0]
                os.system(f"cp {os.path.join(output_result,'rank1_protein.pdb')} {refined_protein_path}")
                os.system(f"cp {os.path.join(output_result,'rank1.sdf')} {refined_sdf_path}")

            except KeyboardInterrupt:
                raise KeyboardInterrupt
            except Exception as e:
                print(f"###### Error: {raw_sdf_path}")
                print(e)
                pass
        
        print(f"Finish {task}")

    def process_items(self, tasks, gpu_id):
        for task in tasks:
            self.process_single_item(task, gpu_id)
    
    def run(self,wd):
        self.wd = wd
        tasks=[f for f in glob.glob(self.wd + "/[0-9A-Z]*_*") if os.path.isdir(f)]
        print(f"Total number of cases: {len(tasks)}")

        task_splits = [tasks[i::len(self.gpu_ids)] for i in range(len(self.gpu_ids))]
        print("tasks_splits",task_splits)

        with concurrent.futures.ProcessPoolExecutor(max_workers=len(self.gpu_ids)) as executor:
            futures = []
            for i in range(len(self.gpu_ids)):
                future = executor.submit(self.process_items, task_splits[i], self.gpu_ids[i])
                futures.append(future) 
            concurrent.futures.wait(futures)

        # for i in range(len(self.gpu_ids)):
        #     self.process_items(task_splits[i])
        print("Finish all cases")

class GninaDockandAssembleRunner:
    def __init__(self, gpu_ids=None, overwrite=False):
        """
        初始化类，设置 GPU IDs 和是否覆盖现有文件。

        :param gpu_ids: 可用的 GPU ID 列表
        :param overwrite: 是否覆盖已有的文件夹
        """
        self.gpu_ids = gpu_ids if gpu_ids is not None else [1, 2, 3, 4, 5, 6, 7]
        self.overwrite = overwrite
        self.name2protein_path = {}

    def get_raw_protein_path(self, pkl_dir, target_name):
        """
        根据目标名称从指定的 .pkl 目录中获取蛋白质文件路径。

        :param pkl_dir: .pkl 文件所在的目录
        :param target_name: 目标的名称
        :return: 蛋白质文件的路径，如果未找到则返回 None
        """
        pkl_files = glob.glob(os.path.join(pkl_dir, "*.pkl"))
        if len(pkl_files) == 0:
            print(f"Error: No .pkl files found in {pkl_dir}")
            return None
        pkl = pkl_files[0]  # 假设只有一个 .pkl 文件
        with open(pkl, "rb") as f:
            try:
                init_items = pickle.load(f)
                for item in init_items:
                    self.name2protein_path[item["pocket_name"]] = item["protein_path"]
            except Exception as e:
                print(f"Error loading .pkl file {pkl}: {e}")
                return None
        return self.name2protein_path.get(target_name)

    def process_single_item(self, task, pkl_dir, gpu_id):
        """
        处理单个任务，包括准备文件、调用 GNINA 进行对接、处理输出结果。

        :param task: 单个任务的目录路径
        :param pkl_dir: .pkl 文件所在的目录路径
        :param gpu_id: 使用的 GPU ID
        """
        target_wd = task
        ligand_files = [f for f in glob.glob(os.path.join(target_wd, "*.sdf"))]

        if not ligand_files:
            print(f"No .sdf files found in {target_wd}. Skipping this task.")
            return

        for raw_sdf_path in tqdm(ligand_files, desc=f"Processing {os.path.basename(target_wd)}"):
            try:
                print("Processing:", raw_sdf_path)
                # 为每个配体创建 mol_path
                mol_path = raw_sdf_path.replace(".sdf", "")
                if self.overwrite and os.path.exists(mol_path):
                    shutil.rmtree(mol_path)
                os.makedirs(mol_path, exist_ok=True)

                # 复制 raw_sdf 到 mol_path
                shutil.copy(raw_sdf_path, mol_path)
                raw_sdf_path = os.path.join(mol_path, os.path.basename(raw_sdf_path))

                # 复制蛋白质文件到 mol_path
                pocket_name = os.path.basename(target_wd).split("_")[0]  # 假设任务目录名格式为 "PocketName_Something"
                source_protein_path = self.get_raw_protein_path(pkl_dir, pocket_name)
                if not source_protein_path:
                    print(f"Error: Protein path not found for target {pocket_name}")
                    continue
                raw_protein_path = os.path.join(mol_path, os.path.basename(source_protein_path).replace("_refined", ""))
                shutil.copy(source_protein_path, raw_protein_path)

                # 如果优化后的蛋白质和配体文件已存在，则跳过
                refined_protein_path = raw_protein_path.replace(".pdb", "_refined.pdb")
                refined_sdf_path = raw_sdf_path.replace(".sdf", "_refined.sdf")
                if os.path.exists(refined_protein_path) and os.path.exists(refined_sdf_path):
                    print(f"Skipping {raw_sdf_path} as output files already exist.")
                    continue

                # 构建 Docker 命令（保持不变）
                docker_command = [
                    "docker", "run", "--rm",
                    f"--gpus", f"device={gpu_id}",
                    "-v", f"{mol_path}:/scr",
                    "gnina/gnina",
                    "gnina",
                    "-r", "/scr/receptor.pdb",
                    "-l", "/scr/crystal_ligand.sdf",
                    "--autobox_ligand", "/scr/crystal_ligand.sdf",
                    "-o", "/scr/docked.sdf",
                    "--flexdist", "3.5",
                    "--flexdist_ligand", "/scr/crystal_ligand.sdf",
                    "--out_flex", "/scr/flex_receptor.pdb",
                    "--full_flex_output",
                    "--num_modes", "2"
                ]

                # 执行 Docker 命令
                docker_command_str = " ".join(docker_command)
                print(f"Running Docker command:\n{docker_command_str}")
                result = subprocess.run(docker_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=False)

                # 检查运行结果
                if result.returncode != 0:
                    print(f"Error running GNINA for {raw_sdf_path}:")
                    print(result.stderr)
                    continue
                else:
                    print(f"Successfully docked {raw_sdf_path}")

                # 处理输出文件
                docked_sdf = os.path.join(mol_path, "docked.sdf")
                flex_receptor_pdb = os.path.join(mol_path, "flex_receptor.pdb")

                if os.path.exists(docked_sdf) and os.path.exists(flex_receptor_pdb):
                    # 重命名文件
                    rank1_sdf = os.path.join(mol_path, "rank1.sdf")
                    rank1_protein_pdb = os.path.join(mol_path, "rank1_protein.pdb")

                    os.rename(docked_sdf, rank1_sdf)
                    os.rename(flex_receptor_pdb, rank1_protein_pdb)

                    # 复制优化后的文件到目标路径
                    shutil.copy(rank1_sdf, refined_sdf_path)
                    shutil.copy(rank1_protein_pdb, refined_protein_path)
                else:
                    print(f"Expected output files not found in {mol_path}")

            except KeyboardInterrupt:
                raise KeyboardInterrupt
            except Exception as e:
                print(f"###### Error processing {raw_sdf_path}")
                print(e)
                pass

        print(f"Finish processing task: {task}")

    def process_items(self, tasks, pkl_dir, gpu_id):
        """
        处理一组任务。

        :param tasks: 任务目录列表
        :param pkl_dir: .pkl 文件所在的目录路径
        :param gpu_id: 使用的 GPU ID
        """
        for task in tasks:
            self.process_single_item(task, pkl_dir, gpu_id)

    def run(self, tasks_dir, pkl_dir):
        """
        启动对接任务，包括任务分配和并行处理。

        :param tasks_dir: 任务目录的根路径
        :param pkl_dir: .pkl 文件所在的目录路径
        """
        self.wd = os.path.expanduser(tasks_dir)
        print(f"Tasks working directory: {self.wd}")
        print(f"PKL files directory: {pkl_dir}")

        # 查找所有符合模式的任务目录
        tasks = [f for f in glob.glob(os.path.join(self.wd, "[0-9A-Z]*_*")) if os.path.isdir(f)]
        print(f"Total number of cases: {len(tasks)}")
        print(f"Tasks list: {tasks}")

        if not tasks:
            print("No tasks found. Please check the tasks directory and ensure it contains the necessary task directories.")
            return

        # 任务分配到各个 GPU
        task_splits = [tasks[i::len(self.gpu_ids)] for i in range(len(self.gpu_ids))]
        print("Tasks splits:", task_splits)

        with concurrent.futures.ProcessPoolExecutor(max_workers=len(self.gpu_ids)) as executor:
            futures = []
            for i in range(len(self.gpu_ids)):
                future = executor.submit(self.process_items, task_splits[i], pkl_dir, self.gpu_ids[i])
                futures.append(future)
            concurrent.futures.wait(futures)

        print("Finish all cases")
        
if __name__ == '__main__':  
    dataset = BfnCrossDockedBackboneCrossDockedUnrelaxed()
    print(dataset)
    runner = RelaxAndAssembleRunner(dataset)
    runner.run()

    runner=GninaDockandAssembleRunner()
    tasks_directory = "/mnt/nfs-ssd/data/itersbdd_bfn/gnina"
    pkl_directory = "/mnt/nfs-ssd/data/itersbdd_bfn/iter_gnina_0"
    runner.run(tasks_directory, pkl_directory)
    
    
    pass
