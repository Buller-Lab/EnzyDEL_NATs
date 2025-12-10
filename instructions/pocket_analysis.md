## Identification of cavities and calculation of their dimensions based on pdb input
Activate conda environment
```bash
conda activate pykvfinder_env
```
### Run the pocket analysis for 1W4T (with adapted probe out and volume cutoffs & 1W4T as reference for superposition)
```bash
python scripts/analyze_pockets.py data/1W4T.pdb --probe_out 8.0 --volume_cutoff 50.0 --reference data/1W4T.pdb
```
Expected Output:
```console
Chain A:
- Protein PDB: results_1W4T/protein_1W4T_A_aligned.pdb
- TOML: results_1W4T/results_1W4T_A.toml
- PDB Output: results_1W4T/output_1W4T_A.pdb
- Histograms: results_1W4T/histograms_1W4T_A.pdf
- Cavities PDB: results_1W4T/cavities_1W4T_A.pdb
Combined results saved: results_1W4T/pocket_results_1W4T.xlsx
```
### Run the pocket analysis for 7QI3 (with adapted probe out and volume cutoffs & 1W4T as reference for superposition)
```bash
python scripts/analyze_pockets.py data/7QI3.pdb --probe_out 8.0 --volume_cutoff 50.0 --reference data/1W4T.pdb
```
