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

git clone https://github.com/Buller-Lab/EnzyDEL_NATs

this should only take a few seconds.

# System Requirements

## Hardware requirements

This code was developed and tested on the following hardware:

CPU: AMD Ryzen Threadripper 3970X 32-Core Processor
Memory: 130 GiB RAM
GPU: 2x NVIDIA GeForce RTX 3090

## Software requirements
We installed the anaconda python distribution from: https://www.anaconda.com/products/individual and followed their download instructions. The specific versions we used for our analysis can be found in requirements.txt


### The followings scripts are provided:
a) pocket_analysis.py 
b) boltz2x_cofolding.py
c) identify_centroid.py



# References

If you utilize this code, please cite:

add reference here
