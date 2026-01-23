import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
import os
import json
from sklearn.cluster import AgglomerativeClustering
import multiprocessing as mp
from tqdm import tqdm
from functools import partial

os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
class ClusterSelect(Select):
    def __init__(self, cluster):
        self.cluster = cluster
    
    def accept_residue(self, residue):
        return residue.get_id()[1] in self.cluster
def domain_split(pdb,pae,output_dir,plddt_thres=50,pae_thres=15):
    #read pdb with biopython, and return the array of plddt,return residue index that plddt > threshold
    parser = PDBParser()
    chain = parser.get_structure('protein',pdb)[0]['A']
    plddt = []
    for residue in chain:
        plddt.append(residue['CA'].get_bfactor())
    plddt = np.array(plddt)
    # Quality filter disabled - process all proteins
    # if plddt.mean() < 69.9 or sum(plddt>=70)/len(plddt) < 0.899:
    #     raise ValueError("Low quality model")
    res_num = np.arange(len(plddt),dtype=int) + 1 #residue index start from 1
    res_num = res_num[plddt >= plddt_thres]

    #parse pae as a josn and return np array of pae
    pae = json.load(open(pae))
    pae = np.array(pae[0]["predicted_aligned_error"])
    pae = np.stack([pae,pae.T]).min(axis=0)

    mask = np.ones_like(pae)
    mask = np.tril(np.triu(mask,-10),10)*35
    pae = pae + mask
    # plt.hist(pae.flatten(),bins=100)
    # plt.savefig(pdb.replace('.pdb',"_hist.png"))
    # plt.close()

    #remove col and rows that plddt < threshold
    pae = pae[plddt >= plddt_thres]
    pae = pae[:,plddt >= plddt_thres]

    # plt.imshow(pae<pae_thres)
    # plt.savefig(pdb.replace('.pdb',".png"))

    #cluster the pae matrix at thres and return the cluster index
    cluster_labels = AgglomerativeClustering(n_clusters=None,metric='precomputed',distance_threshold=pae_thres,linkage='average').fit_predict(pae)
    #for each cluster, get the residue index and write to a list
    clusters = []
    for i in range(max(cluster_labels)+1):
        clusters.append(res_num[cluster_labels == i])
    #sort the clusters by length
    clusters = sorted(clusters,key=lambda x:len(x),reverse=True)
    #for each cluster output pdb files with name of pdbname_clusterindex.pdb
    io = PDBIO()
    output_num = 0
    for i, cluster in enumerate(clusters):
        if len(cluster)<=10:
            break
        #sort the cluster by residue index, remove continous fragment that is shorter than 10
        cluster = sorted(cluster)
        break_point = [0]+[i for i in range(1,len(cluster)) if cluster[i]-cluster[i-1]>1]+[len(cluster)]
        cluster = sum([cluster[break_point[i]:break_point[i+1]] for i in range(len(break_point)-1) if break_point[i+1]-break_point[i]>10],[])
        if len(cluster)<=10:
            continue
        io.set_structure(chain)
        output_pdb = os.path.join(output_dir,os.path.basename(pdb).replace('.pdb',f'_{i}.pdb'))
        io.save(output_pdb, select=ClusterSelect(cluster))
        output_num += 1
    #print(f"output {output_num} domains")
    if output_num == 0:
        raise ValueError("No domain found")
    return output_num
    



def process_pdb(pdbfile, pdb_dir, pae_dir, output_dir):
    try:
        pdb_path = os.path.join(pdb_dir, pdbfile)
        pae_file = pdbfile.strip().replace('-model_v6.pdb', '-predicted_aligned_error_v6.json')
        pae_path = os.path.join(pae_dir, pae_file)
        domain_split(pdb_path, pae_path, output_dir)
        return (pdbfile.strip(), None)
    except Exception as e:
        return (pdbfile.strip(), str(e))

def callback(return_data,log, tbar):
    pdbfile, error = return_data
    tbar.update(1)
    if error is not None:
        log.write(f"{pdbfile}\t{error}\n")
        log.flush()


def main(protein_id_list,total_files,pdb_dir,pae_dir,output_dir,N_CPU,log_path):

    # Read all protein IDs
    f = open(protein_id_list)
    log= open(log_path, 'w')
    tbar = tqdm(total=total_files)
    partial_callback = partial(callback,log = log, tbar=tbar)
    pool = mp.Pool(N_CPU)
    for pdb in f:
        pool.apply_async(process_pdb, args=(pdb.strip(), pdb_dir, pae_dir, output_dir), callback=partial_callback)
    pool.close()
    pool.join()
    f.close()
    log.close()
if __name__ == '__main__':
    # domain_split(
    #     '/data_hdd/home/jiayinjun/Drug-target_complex_cleanup/tmp/AF-A0A4R1UMU8-F1-model_v4.pdb',
    #     '/data_hdd/home/jiayinjun/Drug-target_complex_cleanup/tmp/AF-A0A4R1UMU8-F1-predicted_aligned_error_v4.json',
    #     '/data_hdd/home/jiayinjun/Drug-target_complex_cleanup/tmp/',
    # )
    base_dir = '/home/tanhaichuan/GenPack_for_galaxy/af2db_downloads'
    protein_id_list = os.path.join(base_dir, 'pdb_file_list.txt')
    pdb_dir = os.path.join(base_dir, 'pdb_files')
    pae_dir = os.path.join(base_dir, 'pae_files')
    output_dir = os.path.join(base_dir, 'domains')
    N_CPU = 32  # 根据系统调整
    total_files = 72740
    log_path = os.path.join(base_dir, 'domain_split_errors.txt')
    
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)
    
    main(protein_id_list, total_files, pdb_dir, pae_dir, output_dir, N_CPU, log_path)


    
