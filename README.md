# Pocket Detection part of "Drug The Whole Genome" project

Two methods for pocket detection are implemented in this repository:

## Pocket Detection with Generative model Refinement

We first use [Fpocket](http://fpocket.sourceforge.net/) to detect pockets in size of 10A, while may not be accurate enough. Then we use our SBDD model to generate ligands in the 10A pocket and relax the protein-ligand complex. We use the relaxed ligand to better detect the pocket with size 6A.

- Fpocket code : `./SBDD_AUG/fpocket.py`

- SBDD model code (with the code fork by MolCRAFT): `./SBDD_AUG/sbdd`

- - to use the code, please refer to the [MolCRAFT](https://github.com/AlgoMole/MolCRAFT?tab=readme-ov-file) repository:

- - data for training out SBDD model and the trained checkpoint can be found in google drive: [https://drive.google.com/file/d/1aUGnX_QV3MQKChosdjS0_zXPn7b0YMUP/view?usp=drive_link](https://drive.google.com/file/d/1ddvtLqCAU2qUAxsiO636HWTKVW7Otx10/view?usp=sharing)

- - to use the trained model, please download the checkpoint and put it in the `./SBDD_AUG/sbdd/checkpoints` folder

- Relaxation code: `./SBDD_AUG/relax.py`

- Extraction 6A pocket code with the relaxed ligand: `./SBDD_AUG/extract_pocket.py`

## Pocket Detection with template matching

In this method, we use the [PDBBind](http://www.pdbbind.org.cn/) database as the template database, and use TM-align to align the protein of human. If the TM-score and pocket matching IOU is higher than a threshold, we consider a new pocket is detected.

- TM-align code: `./template_matching/tmalign.py` this code calculates the TM-score between two proteins and save the rotation matrix

- Pocket matching code: `./template_matching/pocket_match.py` this code filters the pocket with the rotation matrix and calculate the IOU

- relax code: `./template_matching/relax.py` this code relax the protein-ligand complex

- remove illegal complex code: `./template_matching/remove_illegal_complex.py` this code removes the illegal complex with the ligand
