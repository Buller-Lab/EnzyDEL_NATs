import mdtraj as md
import numpy as np
import os
import argparse
from sklearn.cluster import KMeans

# Set up argument parser
parser = argparse.ArgumentParser(description="Process .cif files and cluster based on RMSD")
parser.add_argument('--input_folder', required=True, help='Path to input folder containing .cif files')
parser.add_argument('--output_folder', required=True, help='Path to output folder for results')
args = parser.parse_args()

# Get directory paths from arguments
input_directory = args.input_folder
output_directory = args.output_folder

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Step 1: Collect all .cif files in the input directory
cif_files = [f for f in os.listdir(input_directory) if f.endswith('.cif')]

if not cif_files:
    raise FileNotFoundError(f"No .cif files found in {input_directory}.")

# Step 2: Load the first structure as the reference
ref_traj = md.load(os.path.join(input_directory, cif_files[0]))
ref_protein = ref_traj.top.select("protein")  # Select protein atoms for superposition
ref_chain1 = ref_traj.top.select("chainid 1")  # Select chain ID 1 for RMSD

# Step 3: Load and superpose all trajectories to the reference
trajectories = []
for cif_file in cif_files:
    traj = md.load(os.path.join(input_directory, cif_file))
    traj.superpose(ref_traj, atom_indices=ref_protein)  # Superpose by protein atoms
    trajectories.append(traj)

# Step 4: Compute RMSD for chain ID 1
rmsds = []
for traj in trajectories:
    # Compute RMSD for chain ID 1 atoms only, relative to reference chain ID 1
    rmsd = md.rmsd(traj, ref_traj, atom_indices=ref_chain1, ref_atom_indices=ref_chain1)
    rmsds.append(rmsd[0])  # Single frame, so take first element

rmsds = np.array(rmsds).reshape(-1, 1)  # Reshape for clustering

# Step 5: Perform K-means clustering (adjust n_clusters as needed)
n_clusters = min(3, len(cif_files))  # Example: use 3 clusters or fewer if less files
kmeans = KMeans(n_clusters=n_clusters, random_state=42)
kmeans.fit(rmsds)
labels = kmeans.labels_

# Step 6: Find the centroid of the largest cluster
cluster_sizes = np.bincount(labels)
largest_cluster = np.argmax(cluster_sizes)
cluster_indices = np.where(labels == largest_cluster)[0]

# Compute mean RMSD for the largest cluster
cluster_rmsds = rmsds[cluster_indices]
centroid_idx = cluster_indices[np.argmin(np.abs(cluster_rmsds - np.mean(cluster_rmsds)))]

# Step 7: Save the centroid structure (entire system: protein + chain ID 1) as PDB
centroid_traj = trajectories[centroid_idx]
output_file = os.path.join(output_directory, 'centroid.pdb')
centroid_traj.save_pdb(output_file)

print(f"Centroid structure saved as '{output_file}' from file: {cif_files[centroid_idx]}")
