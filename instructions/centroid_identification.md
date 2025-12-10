## How to run the MDTraj & scikit-learn clustering to derive a representative conformation


Clustering of 30 Boltz-2 poses using ligand RMSD and determination of centroid of biggest cluster) 


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
