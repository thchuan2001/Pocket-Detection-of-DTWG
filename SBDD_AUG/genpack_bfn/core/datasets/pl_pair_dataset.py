import os
import pickle
import lmdb
from torch.utils.data import Dataset
from tqdm.auto import tqdm
import sys
from time import time

import torch
from torch_geometric.transforms import Compose

from core.datasets.utils import PDBProtein, parse_sdf_file, ATOM_FAMILIES_ID
from core.datasets.pl_data import ProteinLigandData, torchify_dict

import core.utils.transforms as trans
from multiprocessing import Pool



class DBReader:
    def __init__(self, path) -> None:
        self.path = path
        self.db = None
        self.keys = None

    def _connect_db(self):
        """
            Establish read-only database connection
        """
        assert self.db is None, 'A connection has already been opened.'
        self.db = lmdb.open(
            self.path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=False,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        )
        with self.db.begin() as txn:
            self.keys = list(txn.cursor().iternext(values=False))

    def _close_db(self):
        self.db.close()
        self.db = None
        self.keys = None

    def __del__(self):
        if self.db is not None:
            self._close_db()

    def __len__(self):
        if self.db is None:
            self._connect_db()
        return len(self.keys)

    def __getitem__(self, idx):
        if self.db is None:
            self._connect_db()
        key = self.keys[idx]
        data = pickle.loads(self.db.begin().get(key))
        data = ProteinLigandData(**data)
        data.id = idx
        assert data.protein_pos.size(0) > 0
        return data
    

class PocketLigandPairDataset(Dataset):

    def __init__(self, raw_path, input_data_split_path,output_id_split_path,index_path,dataset_name,split_name, transform=None, version='final'):
        super().__init__()



        self.raw_path = raw_path.rstrip('/')
        self.index_path = os.path.join(index_path, f'{dataset_name}.pkl')
        self.processed_path = os.path.join(index_path, f'{dataset_name}.lmdb')
        self.transform = transform
        self.reader = DBReader(self.processed_path)


        if raw_path.endswith('.pkl'):
            print("load from pkl")
            self.processed_path = os.path.join("/".join(raw_path.split('/')[:-1]), dataset_name + f'.lmdb')
            self.id_split_file=self.processed_path.replace(".lmdb","_id_split.pkl")

            # if processed_path exists and larger than 10M
            if os.path.exists(self.processed_path) and os.path.getsize(self.processed_path) > 10*(1024*1024):
                pass
            else:
                print("process SIU data from pkl")
                self._process_from_pkl()
                print("finish process SIU data from pkl")
            return 


        # if lmdb is smaller than 1M, then remove it
        if os.path.exists(self.processed_path):
            if os.path.getsize(self.processed_path) < 1024 * 1024:
                os.remove(self.processed_path)
                print(f'{self.processed_path} is too small, removed')
                if os.path.exists(self.processed_path + '-lock'):
                    os.remove(self.processed_path + '-lock')
                    print(f'{self.processed_path}-lock is removed')

        new_lmdb=False
        if not os.path.exists(self.processed_path):
            print(f'{self.processed_path} does not exist, begin processing data')
            self._process()
            new_lmdb=True


        if input_data_split_path is not None:
            self.data_split_file = os.path.join(input_data_split_path, split_name + '.pkl')
            self.id_split_file = os.path.join(output_id_split_path, split_name + '.pt')
            if not os.path.exists(self.id_split_file) or new_lmdb:
                if new_lmdb:
                    print("new lmdb, refresh id_split_file")
                else :
                    print(f'{self.id_split_file} does not exist, begin processing id split file')
                self._process_id_split_file()

    def _process_id_split_file(self):
        id_split={
            "train":[],
            "val":[],
            "test":[]
        }

        # split data
        with open(self.data_split_file, 'rb') as f:
            data_split = pickle.load(f)
            
        for i in tqdm(range(len(self))):

            name=self.__getitem__(i).protein_filename
            if name in data_split["train"]:
                id_split["train"].append(i)
            elif name in data_split.get("valid",[]) or name in data_split.get("val",[]) or name in data_split.get("test",[]):
                id_split["test"].append(i)

            # # for crossdocked
            # name=self.__getitem__(i).ligand_filename.replace(".sdf","").split("/")[-1]
            # if name in data_split["train"]:
            #     id_split["train"].append(i)
            # elif name in data_split["valid"]:
            #     id_split["test"].append(i)

        torch.save(id_split, self.id_split_file)
       


    def _process(self):
        db = lmdb.open(
            self.processed_path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=True,
            subdir=False,
            readonly=False,  # Writable
        )
        with open(self.index_path, 'rb') as f:
            index = pickle.load(f)

        num_skipped = 0
        num_data=0
        with db.begin(write=True, buffers=True) as txn:
            for i, (pocket_fn, ligand_fn, *_) in enumerate(tqdm(index)):
                if pocket_fn is None: continue
                try:
                    #data_prefix = '/data/'
                    data_prefix = self.raw_path
                    pocket_dict = PDBProtein(os.path.join(data_prefix, pocket_fn)).to_dict_atom()
                    ligand_dict = parse_sdf_file(os.path.join(data_prefix, ligand_fn))
                    data = ProteinLigandData.from_protein_ligand_dicts(
                        protein_dict=torchify_dict(pocket_dict),
                        ligand_dict=torchify_dict(ligand_dict),
                    )
                    data.protein_filename = pocket_fn
                    data.ligand_filename = ligand_fn
                    data = data.to_dict()  # avoid torch_geometric version issue

                    txn.put(
                        key=f'{num_data:08d}'.encode(),
                        value=pickle.dumps(data)
                    )
                    num_data+=1
                except:
                    num_skipped += 1
                    print('Skipping (%d) %s' % (num_skipped, ligand_fn, ))
                    continue
        db.close()
    
    def __len__(self):
        return len(self.reader)

    def __getitem__(self, idx):
        data = self.reader[idx]
        if self.transform is not None:
            data = self.transform(data)
        return data


def process_item(input_item):
    try:
        if "pocket_name" not in input_item:
            input_item["pocket_name"] = "Unnamed_Molecule"
        pocket_dict = PDBProtein(input_item["pocket_path"]).to_dict_atom()
        ligand_dict = parse_sdf_file(input_item["ref_ligand_path"])
        data = ProteinLigandData.from_protein_ligand_dicts(
            protein_dict=torchify_dict(pocket_dict),
            ligand_dict=torchify_dict(ligand_dict),
        )
        data.protein_filename = input_item["pocket_name"]
        data.ligand_filename = input_item["pocket_name"]
        data = data.to_dict()  # avoid torch_geometric version issue
        return data, None
    except Exception as e:
        raise e
        return None, (input_item["pocket_name"], e)

class PocketLigandPairDatasetByPkl(Dataset):

    def __init__(self,inference_input_pkl, transform=None, version='final'):
        super().__init__()

        self.items=[]
        self.transform = transform
        with open(inference_input_pkl, 'rb') as f:
            inputs = pickle.load(f)
        num_skipped = 0
        num_data=0
        
        # with Pool(processes=30) as pool:
            # results = list(tqdm(pool.imap(process_item, inputs), total=len(inputs)))
        results=[]
        for input_item in tqdm(inputs):
            results.append(process_item(input_item))

        for data, error in results:
            if data is not None:
                self.items.append(data)
                num_data += 1
            else:
                num_skipped += 1
                pocket_name, e = error
                print('Skipping (%d) %s' % (num_skipped, pocket_name))
                print(e)

    
    def __len__(self):
        return len(self.items)

    def __getitem__(self, idx):
        data = self.items[idx]
        if self.transform is not None:
            data = self.transform(data)
        return data



class PocketLigandGeneratedPairDataset(Dataset):

    def __init__(self, raw_path, transform=None, version='4-decompdiff'):
        raise NotImplementedError
        super().__init__()
        self.raw_path = raw_path.rstrip('/')
        self.generated_path = os.path.join('/sharefs/share/sbdd_data/all_results', f'{version}_docked_pose_checked.pt')
        self.processed_path = os.path.join(os.path.dirname(self.raw_path),
                                           os.path.basename(self.raw_path) + f'_processed_{version}.lmdb')
        self.transform = transform
        self.reader = DBReader(self.processed_path)

        if not os.path.exists(self.processed_path):
            print(f'{self.processed_path} does not exist, begin processing data')
            self._process()
        
    def _process(self):
        db = lmdb.open(
            self.processed_path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=True,
            subdir=False,
            readonly=False,  # Writable
        )
        with open(self.generated_path, 'rb') as f:
            results = torch.load(f)

        num_skipped = 0
        with db.begin(write=True, buffers=True) as txn:
            idx = -1
            for i, res in tqdm(enumerate(results), total=len(results)):
                if isinstance(res, dict):
                    res = [res]
                for r in res:
                    idx += 1
                    mol = r["mol"]
                    ligand_fn = r["ligand_filename"]
                    pocket_fn = os.path.join(
                        os.path.dirname(ligand_fn),
                        os.path.basename(ligand_fn)[:-4] + '_pocket10.pdb'
                    )

                    if pocket_fn is None: continue
                    try:
                        data_prefix = self.raw_path
                        pocket_dict = PDBProtein(os.path.join(data_prefix, pocket_fn)).to_dict_atom()
                        ligand_dict = parse_sdf_file(mol)
                        # ligand_dict = parse_sdf_file(os.path.join(data_prefix, ligand_fn))
                        data = ProteinLigandData.from_protein_ligand_dicts(
                            protein_dict=torchify_dict(pocket_dict),
                            ligand_dict=torchify_dict(ligand_dict),
                        )
                        data.protein_filename = pocket_fn
                        data.ligand_filename = ligand_fn
                        data = data.to_dict()  # avoid torch_geometric version issue
                        txn.put(
                            key=str(idx).encode(),
                            value=pickle.dumps(data)
                        )
                    except Exception as e:
                        num_skipped += 1
                        print('Skipping (%d) %s' % (num_skipped, ligand_fn, ), e)
                        continue
        db.close()
    
    def __len__(self):
        return len(self.reader)

    def __getitem__(self, idx):
        data = self.reader[idx]
        if self.transform is not None:
            data = self.transform(data)
        return data


class PocketLigandPairDatasetFromComplex(Dataset):
    def __init__(self, raw_path, transform=None, version='final', radius=10.0):
        super().__init__()
        raise NotImplementedError
        self.raw_path = raw_path.rstrip('/')
        self.index_path = os.path.join(self.raw_path, 'index.pkl')
        base_name = os.path.basename(self.raw_path)
        if 'pocket' in base_name:
            self.processed_path = os.path.join(os.path.dirname(self.raw_path),
                                               os.path.basename(self.raw_path) + f'_processed_{version}.lmdb')
        else:
            self.processed_path = os.path.join(os.path.dirname(self.raw_path),
                                            os.path.basename(self.raw_path) + f'_pocket{radius}_processed_{version}.lmdb')
        self.transform = transform
        self.reader = DBReader(self.processed_path)

        self.radius = radius

        if not os.path.exists(self.processed_path):
            print(f'{self.processed_path} does not exist, begin processing data')
            self._process()

    def _process(self):
        db = lmdb.open(
            self.processed_path,
            map_size=10*(1024*1024*1024),   # 50GB
            create=True,
            subdir=False,
            readonly=False,  # Writable
            max_readers=256,
        )
        with open(self.index_path, 'rb') as f:
            index = pickle.load(f)

        print('Processing data...', 'index', self.index_path, index[0])

        num_skipped = 0
        with db.begin(write=True, buffers=True) as txn:
            for i, (pocket_fn, ligand_fn, *_) in enumerate(tqdm(index)):
                if pocket_fn is None: continue
                try:
                    data_prefix = self.raw_path
                    # clip pocket
                    ligand_dict = parse_sdf_file(os.path.join(data_prefix, ligand_fn))
                    protein = PDBProtein(os.path.join(data_prefix, pocket_fn))
                    selected = protein.query_residues_ligand(ligand_dict, self.radius)
                    pdb_block_pocket = protein.residues_to_pdb_block(selected)
                    pocket_dict = PDBProtein(pdb_block_pocket).to_dict_atom()

                    # pocket_dict = PDBProtein(os.path.join(data_prefix, pocket_fn)).to_dict_atom()
                    # ligand_dict = parse_sdf_file(os.path.join(data_prefix, ligand_fn))
                    
                    data = ProteinLigandData.from_protein_ligand_dicts(
                        protein_dict=torchify_dict(pocket_dict),
                        ligand_dict=torchify_dict(ligand_dict),
                    )
                    data.protein_filename = pocket_fn
                    data.ligand_filename = ligand_fn
                    data = data.to_dict()  # avoid torch_geometric version issue
                    txn.put(
                        key=str(i).encode(),
                        value=pickle.dumps(data)
                    )
                except Exception as e:
                    num_skipped += 1
                    print('Skipping (%d) %s' % (num_skipped, ligand_fn, ), e)
                    with open('skipped.txt', 'a') as f:
                        f.write('Skip %s due to %s\n' % (ligand_fn, e))
                    continue
        db.close()

    def __len__(self):
        return len(self.reader)

    def __getitem__(self, idx):
        data = self.reader[idx]
        if self.transform is not None:
            data = self.transform(data)
        return data
    

class PocketLigandPairDatasetFeaturized(Dataset):
    def __init__(self, raw_path, ligand_atom_mode,inference_input_pkl=None, version='simple'):
        """
        in simple version, only these features are saved for better IO:
            protein_pos, protein_atom_feature, protein_element, 
            ligand_pos, ligand_atom_feature_full, ligand_element
        """

        if version == 'simple':
            self.features_to_save = [
                'protein_pos', 'protein_atom_feature', 'protein_element', 
                'ligand_pos', 'ligand_atom_feature_full', 'ligand_element',
                'protein_filename', 'ligand_filename',
            ]
        else:
            raise NotImplementedError
        
        self.raw_path = raw_path
        self.ligand_atom_mode = ligand_atom_mode
        self.version = version

        if inference_input_pkl is not None:
            self.init_with_pkl(inference_input_pkl)
            return

        #################################################
        # hard code for dataset
        #################################################
        raw_path="/data/pdbbind_2020/combine_set"
        input_data_split_path="/data/rag/data_splits"
        output_id_split_path="/data/pdbbind_2020/bfn_utils/id_split_files"
        index_path="/data/pdbbind_2020/bfn_utils/"
        dataset_name="PDBBind"
        split_name="PDBBind-DUD_E_FLAPP_0.9"

        # raw_path="/data/pdbbind_cb/raw"
        # input_data_split_path="/data/rag/data_splits"
        # output_id_split_path="/data/pdbbind_cb/bfn_utils/id_split_files"
        # index_path="/data/pdbbind_cb/bfn_utils/"
        # dataset_name="PDBBind_cb"
        # split_name="PDBBind-DUD_E_FLAPP_0.9"

        # raw_path="/data/BioLip/raw"
        # input_data_split_path="/data/rag/data_splits"
        # output_id_split_path="/data/BioLip/bfn_utils/id_split_files"
        # index_path="/data/BioLip/bfn_utils/"
        # dataset_name="BioLip"
        # split_name="BioLip-DUD_E_FLAPP_0.9"

        # raw_path="/data/CrossDocked/crossdocked_pocket10"
        # input_data_split_path="/data/rag/data_splits"
        # output_id_split_path="/data/CrossDocked/bfn_utils/id_split_files"
        # index_path="/data/CrossDocked/bfn_utils/"
        # dataset_name="CrossDocked"
        # split_name="CrossDocked-DUD_E_FLAPP_0.6"

        # raw_path='/data/SIU_bfn/final_dic_10a.pkl'
        # index_path='/data/SIU_bfn'
        # input_data_split_path=None
        # output_id_split_path=None
        # id_split_path='/data/SIU_bfn/'
        # dataset_name='SIU_final'
        # split_name=None

        
        ####### Holo SBDD dataset

        ### backbone version
        # raw_path="/data/CrossDocked_backbone/raw"
        # input_data_split_path="/data/CrossDocked_backbone/bfn_utils"
        # output_id_split_path="/data/CrossDocked_backbone/bfn_utils/id_split_files"
        # index_path="/data/CrossDocked_backbone/bfn_utils"
        # dataset_name="CrossDocked_backbone"
        # split_name="CrossDocked_backbone_split_by_name"

        ### full-atom version
        # raw_path="*"
        # input_data_split_path="/data/CrossDocked/bfn_utils"
        # output_id_split_path="/data/CrossDocked/bfn_utils/id_split_files"
        # index_path="/data/CrossDocked/bfn_utils"
        # dataset_name="CrossDocked"
        # split_name="CrossDocked_split_by_name"


        ###### APObind dataset

        ### backbone version
        # raw_path="/data/CrossDocked_backbone/raw"
        # input_data_split_path="/data/CrossDocked_backbone/APObind_bfn_utils"
        # output_id_split_path="/data/CrossDocked_backbone/APObind_bfn_utils/id_split_files"
        # index_path="/data/CrossDocked_backbone/APObind_bfn_utils"
        # dataset_name="CrossDocked_backbone_APObind"
        # split_name="CrossDocked_backbone_APObind_split_by_name"

        ### full-atom version
        # raw_path = "*"
        # input_data_split_path = "/data/CrossDocked/APObind_bfn_utils"
        # output_id_split_path = "/data/CrossDocked/APObind_bfn_utils/id_split_files"
        # index_path = "/data/CrossDocked/APObind_bfn_utils"
        # dataset_name = "CrossDocked_APObind"
        # split_name = "CrossDocked_APObind_split_by_name"



        print("#"*30)
        print("Hard code for dataset: please press Y to continue")
        print("#"*30)
        print(f"\traw_path: {raw_path}")
        print(f"\tinput_data_split_path: {input_data_split_path}")
        print(f"\toutput_id_split_path: {output_id_split_path}")
        print(f"\tindex_path: {index_path}")
        print(f"\tdataset_name: {dataset_name}")
        print(f"\tsplit_name: {split_name}")
        print("#"*30)
        assert input("Continue? (Y/N): ").lower() == "y" , "stop"

        #################################################
        # end of hard code for dataset
        #################################################

        self.raw_path=raw_path
        self.input_data_split_path=input_data_split_path
        self.output_id_split_path=output_id_split_path
        self.index_path=index_path
        self.dataset_name=dataset_name
        if split_name is None:
            split_name=self.dataset_name
        self.split_name=split_name 

        # self.transformed_path = os.path.join(
        #     os.path.dirname(self.raw_path), os.path.basename(self.raw_path) + 
        #     f'_{ligand_atom_mode}_transformed_{version}.pt'
        # )
        self.transformed_path = os.path.join(self.index_path,split_name + f'_{ligand_atom_mode}_transformed_{version}.pt')
        if not os.path.exists(self.transformed_path):
            print(f'{self.transformed_path} does not exist, begin transforming data')
            self._transform()
        else:
            print(f'reading data from {self.transformed_path}...')
            tic = time()
            tr_data = torch.load(self.transformed_path)
            toc = time()
            print(f'{toc - tic} elapsed')
            # debug
            self.train_data, self.val_data, self.test_data = tr_data['train'], tr_data['val'][:100], tr_data['test'][:3]
            self.protein_atom_feature_dim = tr_data['protein_atom_feature_dim']
            self.ligand_atom_feature_dim = tr_data['ligand_atom_feature_dim'] 
        
    def _transform(self):
        raw_dataset = PocketLigandPairDataset(self.raw_path,
                                              self.input_data_split_path,
                                                self.output_id_split_path,
                                                self.index_path,
                                                self.dataset_name,
                                                self.split_name,
                                               None, 'final')
        
        if self.raw_path.endswith('.pkl'):
            processed_path = os.path.join("/".join(self.raw_path.split('/')[:-1]), self.dataset_name + f'.lmdb')
            split_path=processed_path.replace(".lmdb","_id_split.pkl")
        else :
            split_path = os.path.join(self.output_id_split_path, self.split_name + '.pt')

        if split_path.endswith('.pkl'):
            split = pickle.load(open(split_path, 'rb'))
        elif split_path.endswith('.pt'):
            split = torch.load(split_path) 
        else:
            raise ValueError(f'Unknown split file type: {split}')
        
        if 'val' not in split or len(split['val']) == 0:
            split['val'] = split['test']
        train_ids, val_ids, test_ids = split['train'], split['val'], split['test']
        print(f'train_size: {len(train_ids)}, val_size: {len(val_ids)}, test_size: {len(test_ids)}')

        protein_featurizer = trans.FeaturizeProteinAtom()
        ligand_featurizer = trans.FeaturizeLigandAtom(self.ligand_atom_mode)
        transform_list = [
            protein_featurizer,
            ligand_featurizer,
            # trans.FeaturizeLigandBond(),
        ]
        transform = Compose(transform_list)
        self.protein_atom_feature_dim = protein_featurizer.feature_dim
        self.ligand_atom_feature_dim = ligand_featurizer.feature_dim

        def _transform_subset(ids):
            data_list = []

            for idx in tqdm(ids):
                data = raw_dataset[idx]
                data = transform(data)
                tr_data = {}
                for k in self.features_to_save:
                    tr_data[k] = getattr(data, k)
                tr_data['id'] = idx
                tr_data = ProteinLigandData(**tr_data)
                data_list.append(tr_data)
            return data_list

        self.train_data = _transform_subset(train_ids)
        print(f'train_size: {len(self.train_data)}, {sys.getsizeof(self.train_data)}')
        self.val_data = _transform_subset(val_ids)
        print(f'val_size: {len(self.val_data)}, {sys.getsizeof(self.val_data)}')
        self.test_data = _transform_subset(test_ids)
        print(f'test_size: {len(self.test_data)}, {sys.getsizeof(self.test_data)}')
        torch.save({
            'train': self.train_data, 'val': self.val_data, 'test': self.test_data, 
            'protein_atom_feature_dim': self.protein_atom_feature_dim,
            'ligand_atom_feature_dim': self.ligand_atom_feature_dim,
        }, self.transformed_path)

    def init_with_pkl(self,inference_input_pkl):
        raw_dataset = PocketLigandPairDatasetByPkl(inference_input_pkl, None, 'final')

        protein_featurizer = trans.FeaturizeProteinAtom()
        ligand_featurizer = trans.FeaturizeLigandAtom(self.ligand_atom_mode)
        transform_list = [
            protein_featurizer,
            ligand_featurizer,
            # trans.FeaturizeLigandBond(),
        ]
        transform = Compose(transform_list)
        self.protein_atom_feature_dim = protein_featurizer.feature_dim
        self.ligand_atom_feature_dim = ligand_featurizer.feature_dim

        def _transform_subset(ids):
            data_list = []

            for idx in tqdm(ids):
                data = raw_dataset[idx]
                data = ProteinLigandData(**data)

                data = transform(data)
                tr_data = {}
                for k in self.features_to_save:
                    tr_data[k] = getattr(data, k)
                tr_data['id'] = idx
                tr_data = ProteinLigandData(**tr_data)
                data_list.append(tr_data)
            return data_list

        self.train_data = []
        print(f'train_size: {len(self.train_data)}, {sys.getsizeof(self.train_data)}')
        self.val_data = []
        print(f'val_size: {len(self.val_data)}, {sys.getsizeof(self.val_data)}')
        self.test_data = _transform_subset(range(len(raw_dataset)))
        print(f'test_size: {len(self.test_data)}, {sys.getsizeof(self.test_data)}')


if __name__ == '__main__':
    # original dataset
    dataset = PocketLigandPairDataset('./data/crossdocked_v1.1_rmsd1.0_pocket10')
    print(len(dataset), dataset[0])

    ############################################################

    # test DecompDiffDataset
    # dataset = PocketLigandGeneratedPairDataset('/sharefs/share/sbdd_data/crossdocked_pocket10')闹闹。
    # print(len(dataset), dataset[0])

    ############################################################

    # test featurized dataset (GPU accelerated)
    # path = '/sharefs/share/sbdd_data/crossdocked_v1.1_rmsd1.0_pocket10'
    # ligand_atom_mode = 'add_aromatic'

    # dataset = PocketLigandPairDatasetFeaturized(path, ligand_atom_mode)
    # train_data, test_data = dataset.train_data, dataset.test_data
    # print(f'train_size: {len(train_data)}, {sys.getsizeof(train_data)}')
    # print(f'test_size: {len(test_data)}, {sys.getsizeof(test_data)}')
    # print(test_data[0], sys.getsizeof(test_data[0]))

    ############################################################

    # test featurization
    # find all atom types
    # atom_types = {(1, False): 0}

    # dataset = PocketLigandPairDataset(path, transform)
    # for i in tqdm(range(len(dataset))):
    #     data = dataset[i]
    #     element_list = data.ligand_element
    #     hybridization_list = data.ligand_hybridization
    #     aromatic_list = [v[trans.AROMATIC_FEAT_MAP_IDX] for v in data.ligand_atom_feature]

    #     types = [(e, a) for e, h, a in zip(element_list, hybridization_list, aromatic_list)]
    #     for t in types:
    #         t = (t[0].item(), bool(t[1].item()))
    #         if t not in atom_types:
    #             atom_types[t] = 0
    #         atom_types[t] += 1

    # idx = 0
    # for k in sorted(atom_types.keys()):
    #     print(f'{k}: {idx}, # {atom_types[k]}')
    #     idx += 1

    ############################################################
    
    # count atom types
    # type_counter, aromatic_counter, full_counter = {}, {}, {}
    # for i, data in enumerate(tqdm(dataset)):
    #     element_list = data.ligand_element
    #     hybridization_list = data.ligand_hybridization
    #     aromatic_list = [v[trans.AROMATIC_FEAT_MAP_IDX] for v in data.ligand_atom_feature]
    #     flag = False

    #     for atom_type in element_list:
    #         atom_type = int(atom_type.item())
    #         if atom_type not in type_counter:
    #             type_counter[atom_type] = 0
    #         type_counter[atom_type] += 1

    #     for e, a in zip(element_list, aromatic_list):
    #         e = int(e.item())
    #         a = bool(a.item())
    #         key = (e, a)
    #         if key not in aromatic_counter:
    #             aromatic_counter[key] = 0
    #         aromatic_counter[key] += 1

    #         if key not in trans.MAP_ATOM_TYPE_AROMATIC_TO_INDEX:
    #             flag = True

    #     for e, h, a in zip(element_list, hybridization_list, aromatic_list):
    #         e = int(e.item())
    #         a = bool(a.item())
    #         key = (e, h, a)
    #         if key not in full_counter:
    #             full_counter[key] = 0
    #         full_counter[key] += 1
        
    # print('type_counter', type_counter)
    # print('aromatic_counter', aromatic_counter)
    # print('full_counter', full_counter)

