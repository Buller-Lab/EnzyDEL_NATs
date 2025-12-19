# Cluster-Model CREST Workflow for Enzyme–Ligand Binding Energies

This folder contains scripts, inputs, instructions, and results to reproduce the cluster-model-based
Distribution-wheighted binding energies calculation workflow. 

---

## Overview and Workflow

The protocol consists of the following steps:

1. Cluster Model (CM) construction from enzyme-binder complexes
2. Semiempirical (GFN2-xtb) geometry optimization of the CM
3. Conformational sampling of the bound binder using CREST in GFNFF level of theory
4. Re-ranking of conformers using higher-accuracy semiempirical methods
5. Boltzmann-weighted binding energy calculation

Schematic overview:

Enzyme–Ligand Structure  
→ Cluster Model (6 Å cutoff)  
→ GFN2-xTB Minimization (fixed backbone)  
→ CREST (GFNFF, -nci, iMTD-GC)  
→ GFN2-xTB / g-xTB Single Points  
→ Boltzmann-Weighted Binding Energy  

---

## Cluster Model Construction

Cluster models are generated starting from enzyme–binder structures.

Protocol:

- All residues with at least one atom within 6 Å of any ligand atom are selected.
- Each residue is truncated beyond the Cα (CA) atom.
- Terminal O, N, C, and H atoms beyond the CA are removed.
- Additional hydrogen atoms are added to the CA to saturate valence.
- CA atoms and the newly added hydrogens are fixed in all subsequent calculations
  to prevent unrealistic backbone distortions.

Script:

- `Pymol_CM_creation.py`  
  Automates truncation, addition of missing H atoms, and preparation of cluster models using PyMOL.

---

## Geometry Optimization

The full cluster model (protein + ligand) is optimized at the GFN2-xTB level of theory.

Details:

- Implicit solvation: ALPB model with water dielectric constant
- Optimizer: Cartesian Fast Inertial Relaxation Engine (FIRE)
- CA atoms and added hydrogens remain fixed
- Final GFN2-xTB atomic charges are used for subsequent GFNFF calculations

References:

- GFN2-xTB: https://doi.org/10.1021/acs.jctc.8b01176  
- ALPB: https://doi.org/10.1021/acs.jctc.1c00471  
- FIRE optimizer: https://doi.org/10.1103/PhysRevLett.97.170201  

---

## Conformational Sampling with CREST

Binding conformers are generated using CREST at the GFNFF level of theory.

Key settings:

- CREST using the iMTD-GC algorithm (default)
- `-nci` keyword to enhance noncovalent interaction sampling
- The protein CM, CA atoms, and added hydrogens remain fixed
- Only selected atoms (ligand, decoy atoms, or DNA-tag atoms up to the second base)
  are allowed to move and participate in the internal metadynamics


For DNA-tag-containing systems:

- DNA is truncated once it no longer interacts with the enzyme (typically the 6th base
  is excluded).
- The terminal phosphate group is removed and the terminal O is protonated and fixed.
- This reduces computational cost and improves charge stability.

References:

- CREST: https://doi.org/10.1063/5.0197592  
- GFNFF: https://doi.org/10.1002/anie.202004239  
- iMTD-GC: https://doi.org/10.1039/C9CP06869D  

---

## Conformer Re-ranking

The final CREST conformers are re-ranked using single-point energy calculations at:

- GFN2-xTB
or
- g-xTB (https://doi.org/10.26434/chemrxiv-2025-bjxvt)

Script:

- `rerank_xtb.py`  
  Performs single-point calculations and energy-based sorting of CREST conformers.
  This is important for the Binding energy calculation, because E_tot will not be re-computed.

---

## Binding Energy Calculation

Binding energies are computed using a Boltzmann-weighted ensemble approach.

Binding energy Definition:

E_bind = E_tot − (E_prot + E_sub)

Where:

- E_tot is the electronic energy of the full CM (protein + binder)
- E_prot is the electronic energy of the CM with the binder removed
- E_sub is the electronic energy of the binder alone

Boltzmann weighting:

- Conformer probabilities are computed using Boltzmann distribution, using E_tot
- Only conformers accounting for 99.9999% of the total population are included
- The final binding energy is a population-weighted average over this ensemble

Script:

- `binding_energies.py`  
  Computes Boltzmann weights and ensemble-averaged binding energies.

---

## Software Requirements

- xTB (GFN2-xTB, GFNFF, g-xTB) (https://xtb-docs.readthedocs.io)
- CREST 3 (https://crest-lab.github.io)
- PyMOL
- Python ≥ 3.8
  - numpy
  - scipy

---

## How to use

- Create the CM. If you plan to use the included pymol script:
1. Select the aminoacids that will be forming the CM. Recomendation: Select the binder and in (sele) `Modify > Arround > residues arround 6Å` or more.
2. In Pymol: `run Pymol_CM_creation.py` and `make_cluster_model sele`.

This will create the proteic part of the cluster model and also print the atoms to fix in further calculations.
Add the binding molecule to the same xyz structure and save it. Make sure that the atom numbers are not shifted (it is recomended to put the binder after the proteic part to avoid atom ID shifting)

- Set threads and memory that xtb will use. This is important or the xtb calculation will not finish (large number of atoms). Ajust for your machine, this are just recomended values:

```
ulimit -s unlimited
export OMP_STACKSIZE=6G
export OMP_NUM_THREADS=16,1
```

Minimize the structure at GFN2-xtb level using xtb.
```
xtb CM_binder.xyz --alpb water --charge -4 --input XTB_constraints.inp --opt
```

Remember to ajust the charge value (`--charge -4` in this example and the list of atoms to fix in the XTB_constraints.inp (example and constraints are available in the xtb_constraint_files folder)

This command will generate a `xtbolt.xyz` file and a `charge` file with and others. I recomend to rename the `xtbopt.xyz` to something that can be easily recognizable (`CM_binder_GFN2opt.xyz`), and `charge` must be placed, without changing the name, in the same place the CREST calculation will take place.

- Compute the conformers using CREST like this.
```
crest crest_input.toml --nci --noopt —cinp CREST_constraints.inp
```

CREST toml input and constraint files are available in the `CREST_inputs_and_constraints` folder. Remember to ajust the atoms fixed, the atoms selected for the metadynamics (binder/part of the binder), charge, and number of threads.

This command will use CREST to create the conformers structures (`crest_conformers.xyz`) for your system that will be later used for the binding energy calculation.

- Compute the energy and rank the conformers:
```
python rerank_xtb.py gfn2 crest_conformers.xyz crest_conformers_gfn2.xyz -4
```  

Ajust the method (`gfn2` or `gxtb`) and charge to suit your system.

- Compute the unbound energies and compute the weighted binding energies.
Examples:
```
python binding_energies.py crest_conformers_GFN2.xyz --engine gfn2 --q_sub 0 --q_prot 3 --sub_sel :60
python binding_energies.py crest_conformers_gxtb.xyz --engine gxtb --q_sub 0 --q_prot 3 --sub_sel 40:60
python binding_energies.py crest_conformers_GFN2.xyz --engine gfn2 --q_sub -7 --q_prot 3 --sub_sel 682:
```

Ajust the selection of method, atoms from the binder and charges of both protein (--q_prot) and binder (--q_sub) accordingly. The total energy of the system will not be re-computed, so it is important that the ranked conformers using the same method is used as input.
