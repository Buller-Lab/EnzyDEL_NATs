import mdtraj as md
import numpy as np
import os
import argparse
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from scipy.optimize import linear_sum_assignment
import matplotlib.pyplot as plt
import seaborn as sns

# ----------------------------- Args -----------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--input_folder', required=True, help="Folder containing .cif files")
parser.add_argument('--output_folder', required=True, help="Folder to save results")
parser.add_argument('--perform_sensitivity_analysis', action='store_true',
                    help="If provided, perform sensitivity analysis on subset sizes")
args = parser.parse_args()

input_dir = args.input_folder
output_dir = args.output_folder
os.makedirs(output_dir, exist_ok=True)

# ----------------------------- Load & align -----------------------------
cif_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.cif')])
if not cif_files:
    raise FileNotFoundError("No .cif files found in the input folder.")
print(f"Found {len(cif_files)} structures.")

ref = md.load(os.path.join(input_dir, cif_files[0]))
ref_protein = ref.top.select("protein")
ref_chain1 = ref.top.select("chainid 1 or chainid 2")

trajectories = [
    md.load(os.path.join(input_dir, f)).superpose(ref, atom_indices=ref_protein)
    for f in cif_files
]

# ----------------------------- RMSD -----------------------------
rmsds = np.array([
    md.rmsd(t, ref, atom_indices=ref_chain1)[0] for t in trajectories
]).reshape(-1, 1)

# ----------------------------- Silhouette -----------------------------
max_k = min(10, len(cif_files) - 1)
sil_data = []
for k in range(2, max_k + 1):
    labels = KMeans(n_clusters=k, random_state=42, n_init=10).fit_predict(rmsds)
    score = silhouette_score(rmsds, labels)
    sil_data.append({"n_clusters": k, "silhouette_score": score})

sil_df = pd.DataFrame(sil_data)
sil_df.to_csv(os.path.join(output_dir, "silhouette_scores.csv"), index=False)

best_k = int(sil_df.loc[sil_df["silhouette_score"].idxmax(), "n_clusters"])
print(f"Best k = {best_k}")

# ----------------------------- Final clustering -----------------------------
km_full = KMeans(n_clusters=best_k, random_state=42, n_init=10).fit(rmsds)
full_labels = km_full.labels_
full_centroids = km_full.cluster_centers_.flatten()

# Save centroid structure (closest to mean in the largest cluster)
largest_cluster = np.argmax(np.bincount(full_labels))
cluster_rmsds = rmsds[full_labels == largest_cluster].flatten()
best_idx = np.where(full_labels == largest_cluster)[0][
    np.argmin(np.abs(cluster_rmsds - cluster_rmsds.mean()))
]
trajectories[best_idx].save_pdb(os.path.join(output_dir, "centroid.pdb"))

pd.DataFrame({
    "filename": cif_files,
    "RMSD_to_ref": rmsds.flatten(),
    "cluster": full_labels
}).to_csv(os.path.join(output_dir, "cluster_assignments.csv"), index=False)

# ----------------------------- Sensitivity Analysis (Optional) -----------------------------
if args.perform_sensitivity_analysis:
    print("Performing sensitivity analysis...")

    def match_shift(c1, c2):
        cost = np.abs(c1.reshape(-1, 1) - c2.reshape(1, -1))
        ri, ci = linear_sum_assignment(cost)
        return cost[ri, ci].mean()

    max_n = min(30, len(cif_files))
    records = []
    for n in range(1, max_n + 1):
        rec = {"num_models": n}
        if n < best_k:
            rec.update({"ari": np.nan, "nmi": np.nan, "centroid_shift": np.nan})
        else:
            sub_rmsds = rmsds[:n]
            km = KMeans(n_clusters=best_k, random_state=42, n_init=10).fit(sub_rmsds)
            rec["ari"] = adjusted_rand_score(full_labels[:n], km.labels_)
            rec["nmi"] = normalized_mutual_info_score(full_labels[:n], km.labels_)
            rec["centroid_shift"] = match_shift(km.cluster_centers_.flatten(), full_centroids)
        records.append(rec)

    sens_df = pd.DataFrame(records)
    sens_df.to_csv(os.path.join(output_dir, "sensitivity_analysis.csv"), index=False)

    # ----------------------------- Plotting -----------------------------
    sns.set_style("whitegrid")
    plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "axes.titlesize": 15})
    fig = plt.figure(figsize=(15, 4.8))
    gs = fig.add_gridspec(1, 3, wspace=0.33)

    # 1. Silhouette
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(sil_df["n_clusters"], sil_df["silhouette_score"],
             marker='o', color='#2c7bb6', linewidth=3, markersize=8)
    ax1.axvline(best_k, color='red', linestyle='--', linewidth=2, label=f'k = {best_k}')
    ax1.set_ylim(0, 1)
    ax1.set_xlabel("Number of clusters (k)")
    ax1.set_ylabel("Silhouette score")
    ax1.set_title("Silhouette Analysis")
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # 2. Agreement scores
    ax2 = fig.add_subplot(gs[0, 1])
    valid = sens_df.dropna()
    ax2.plot(valid["num_models"], valid["ari"], 's-', color='#d7191c', linewidth=3, label="ARI")
    ax2.plot(valid["num_models"], valid["nmi"], '^-', color='#fdae61', linewidth=3, label="NMI")
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("Number of models used")
    ax2.set_ylabel("Agreement score")
    ax2.set_title("Consistency with Full Clustering")
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    # 3. Centroid shift
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.plot(valid["num_models"], valid["centroid_shift"],
             'D-', color='#1a9850', linewidth=3, markersize=7)
    ax3.set_ylim(0, 2)
    ax3.set_xlabel("Number of models used")
    ax3.set_ylabel("Mean centroid shift (nm)")
    ax3.set_title("Centroid Stability")
    ax3.grid(True, alpha=0.3)

    # Save figures
    plt.savefig(os.path.join(output_dir, "clustering_summary.png"), dpi=300, bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, "clustering_summary.pdf"), bbox_inches='tight')
    plt.close()

    print("\nSensitivity analysis completed.")
    print("Figure with fixed axes saved as:")
    print(" clustering_summary.png")
    print(" clustering_summary.pdf")
else:
    print("\nSensitivity analysis skipped (add --perform_sensitivity_analysis to enable).")
    print("Basic clustering results saved.")

print("\nAll done!")
