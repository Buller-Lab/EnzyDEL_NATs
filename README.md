## This is the GitHub Repository for the publication: 

# Engineering an enzyme cascade for a broad substrate scope facilitates DNA-encoded library technology

In this repository we provide the code that was used for binding pocket analysis using pyKVFinder [1] as well as Boltz-2 [2,3] based cofolding of different N-acetyltransferases with their DNA-tagged substrates (30 diffussion model replicates). Subsequently, a ligand RMSD-based clustering script based on MDTraj [4] and scikit-learn [5] was utilized to identify the largest cluster and select it's centroid as representative conformation.

[1] Guerra, João Victor da Silva, et al. "pyKVFinder: an efficient and integrable Python package for biomolecular cavity detection and characterization in data science." BMC bioinformatics 22.1 (2021): 607.

[2] Wohlwend, Jeremy, et al. "Boltz-1 democratizing biomolecular interaction modeling." BioRxiv (2025): 2024-11.

[3] Passaro, Saro, et al. "Boltz-2: Towards accurate and efficient binding affinity prediction." BioRxiv (2025): 2025-06.

[4] McGibbon, Robert T., et al. "MDTraj: a modern open library for the analysis of molecular dynamics trajectories." Biophysical journal 109.8 (2015): 1528-1532.

[5] Pedregosa, Fabian, et al. "Scikit-learn: Machine learning in Python." the Journal of machine Learning research 12 (2011): 2825-2830.

# Installation requirements

this repository can be downloaded to your local machine via the command:
```bash
git clone https://github.com/Buller-Lab/EnzyDEL_NATs
```
then navigate into the cloned repository with:
```bash
cd EnzyDEL_NATs
```

this should only take a few seconds.

# System Requirements

## Hardware requirements

This code was developed and tested on the following hardware:

- CPU: AMD Ryzen Threadripper 3970X 32-Core Processor
- Memory: 130 GiB RAM
- GPU: 2x NVIDIA GeForce RTX 3090

## Software requirements
To create conda environments with necessary dependencies, run:
```bash
conda env create --file pykvfinder_env.yml
```
```bash
conda create -n cofolding_env
git clone https://github.com/jwohlwend/boltz.git
cd boltz; pip install -e .[cuda]
cd ..
```
```bash
conda env create --file clustering_env.yml
```
# Instructions for use
## The followings scripts are provided:
- pocket_analysis.py (identification of cavities and calculation of their dimensions based on pdb input)
- boltz2x_cofolding.py (Boltz-2 co-folding with 30 diffusion models per input yml)
- identify_centroid.py (Clustering of 30 Boltz-2 poses using ligand RMSD and determination of centroid of biggest cluster) 
## How to run the pocket analysis
Activate conda environment
```bash
conda activate pykvfinder_env
```
Run the pocket analysis (with default probe out and volume cutoffs)
```bash
python analyze_pockets.py 1W4T.pdb --probe_out 8.0 --volume_cutoff 50.0
```
```bash
python analyze_pockets.py 7QI3.pdb --probe_out 8.0 --volume_cutoff 50.0
```
## How to run the Boltz-2 cofolding (with 30 diffusions)
Activate conda environment
```bash
conda activate cofolding_env
```
Run the cofolding (with input and output folder specified)
```bash
python boltz2x_cofolding.py --input_folder cofolding_inputs --output_folder cofolding_outputs
```

## How to run the MDTraj & scikit-learn clustering to derive a representative conformation
Activate conda environment
```bash
conda activate clustering_env
```
Run the clustering (with input and output folder specified; only one example is shown and path needs to be adapted for the other variants)
```bash
python identify_centroid.py --input_folder cofolding_outputs/boltz_results_05PaAT_chimera_Substrate/predictions/05PaAT_chimera_Substrate --output_folder cofolding_outputs/centroid_05PaAT_chimera_Substrate
```

# References

If you utilize this code, please cite:

add reference here
