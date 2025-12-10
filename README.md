## This is the GitHub Repository for the publication: 

# A tailored enzyme cascade facilitates DNA-encoded library technology and gives access to a broad substrate scope

In this repository we provide the code that was used for binding pocket analysis using pyKVFinder [1] as well as Boltz-2 [2,3] based cofolding of different N-acetyltransferases with their DNA-tagged substrates (30 diffussion model replicates). Subsequently, a ligand RMSD-based clustering script based on MDTraj [4] and scikit-learn [5] was utilized to identify the largest cluster and select it's centroid as representative conformation.

[1] Guerra, João Victor da Silva, et al. "pyKVFinder: an efficient and integrable Python package for biomolecular cavity detection and characterization in data science." BMC bioinformatics 22.1 (2021): 607.

[2] Wohlwend, Jeremy, et al. "Boltz-1 democratizing biomolecular interaction modeling." BioRxiv (2025): 2024-11.

[3] Passaro, Saro, et al. "Boltz-2: Towards accurate and efficient binding affinity prediction." BioRxiv (2025): 2025-06.

[4] McGibbon, Robert T., et al. "MDTraj: a modern open library for the analysis of molecular dynamics trajectories." Biophysical journal 109.8 (2015): 1528-1532.

[5] Pedregosa, Fabian, et al. "Scikit-learn: Machine learning in Python." the Journal of machine Learning research 12 (2011): 2825-2830.

# Installation

We recommend to run this code on UNIX based systems such as Ubuntu. This repository can be downloaded to your local machine via the command:
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
conda create -n cofolding_env python=3.11 -y
conda activate cofolding_env
git clone https://github.com/jwohlwend/boltz.git
cd boltz
pip install -e ."[cuda]"
cd ..
```
```bash
conda env create --file clustering_env.yml
```
# Instructions for use
## The followings scripts are provided:
- pocket_analysis.py (see [Instructions Pocket Analysis](instructions/pocket_analysis.md); identification of cavities and calculation of their dimensions based on pdb input)
- boltz2x_cofolding.py (Boltz-2 co-folding with 30 diffusion models per input yml)
- identify_centroid.py (Clustering of 30 Boltz-2 poses using ligand RMSD and determination of centroid of biggest cluster) 

## How to run the Boltz-2 cofolding (with 30 diffusions)
Activate conda environment
```bash
conda activate cofolding_env
```
Run the cofolding (with input and output folder specified)
```bash
python scripts/boltz2x_cofolding.py --input_folder cofolding_inputs
 
```
Expected Output: Boltz-2 results folder for every yaml input

## How to run the MDTraj & scikit-learn clustering to derive a representative conformation
Activate conda environment
```bash
conda activate clustering_env
```
### Quick run (clustering only)
 ```bash
 python scripts/identify_centroids.py \
   --input_folder  path/to/your/predictions_folder \
   --output_folder path/to/output_folder
 ```
#
### Full run with sensitivity analysis (check for stability depending on number of Boltz models used)
 ```bash
 python scripts/identify_centroids.py \
   --input_folder  path/to/your/predictions_folder \
   --output_folder path/to/output_folder \
   --perform_sensitivity_analysis
 ```
#
### Example
 ```bash
 python scripts/identify_centroids.py \
   --input_folder boltz_results_42GmAT_Substrate/predictions/42GmAT_Substrate \
   --output_folder centroid_42GmAT_Substrate \
   --perform_sensitivity_analysis
 ```
#
### Output files (always created)
 - `centroid.pdb`                  → Most representative structure
 - `cluster_assignments.csv`       → Filename, RMSD to reference, cluster label
 - `silhouette_scores.csv`         → Justification of chosen number of clusters
#
### Additional files (only when `--perform_sensitivity_analysis` is used)
 - `sensitivity_analysis.csv`      → ARI, NMI, and centroid shift vs. number of models
 - `clustering_summary.png` / `.pdf` → 3-panel figure (silhouette + consistency + centroid stability)
#

# References

If you utilize this code, please cite:
Daniela Schaub, Alice Lessing et al. A tailored enzyme cascade facilitates DNA-encoded library technology and gives access to a broad substrate scope, 23 September 2025, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-7598475/v1]
