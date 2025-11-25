import os
import argparse
import pyKVFinder
import pandas as pd
import mdtraj as md
from pymol import cmd

def align_structures_pymol(mobile_pdb, reference_pdb):
    """
    Align mobile structure to reference using PyMOL's super command.
    
    Args:
        mobile_pdb: Path to mobile PDB file
        reference_pdb: Path to reference PDB file
    
    Returns:
        Path to aligned mobile PDB file
    """
    # Initialize PyMOL in quiet mode
    cmd.reinitialize()
    
    # Load structures
    cmd.load(reference_pdb, "ref")
    cmd.load(mobile_pdb, "mobile")
    
    # Run PyMOL's super alignment
    cmd.super("mobile", "ref")
    
    # Save aligned mobile structure
    aligned_pdb = mobile_pdb.replace('.pdb', '_aligned.pdb')
    cmd.save(aligned_pdb, "mobile")
    
    # Clean up
    cmd.delete("all")
    
    return aligned_pdb

def process_pdb_file(filepath, probe_out, volume_cutoff, reference_pdb=None):
    """Process PDB file with pyKVFinder using protein only, split by chains."""
    filename = os.path.splitext(os.path.basename(filepath))[0]
    output_folder = f'results_{filename}'
    os.makedirs(output_folder, exist_ok=True)
    
    # Load trajectory and select protein atoms only (excluding water)
    traj = md.load(filepath)
    protein_indices = traj.topology.select("protein")
    protein_traj = traj.atom_slice(protein_indices)
    
    # Get unique chains
    chains = set([traj.topology.atom(i).residue.chain.index for i in protein_indices])
    chains = sorted(list(chains))
    
    # Validate reference if provided
    if reference_pdb and not os.path.exists(reference_pdb):
        print(f"Error: Reference PDB file '{reference_pdb}' does not exist")
        return
    
    all_results = []
    
    # Process each chain
    for chain_idx in chains:
        # Select atoms for this chain
        chain_selection = f"protein and chainid {chain_idx}"
        chain_indices = protein_traj.topology.select(chain_selection)
        chain_traj = protein_traj.atom_slice(chain_indices)
        
        # Get chain letter
        chain_letter = traj.topology.atom(protein_indices[0]).residue.chain.chain_id
        for i in protein_indices:
            if traj.topology.atom(i).residue.chain.index == chain_idx:
                chain_letter = traj.topology.atom(i).residue.chain.chain_id
                break
        
        # Save chain-specific protein PDB
        chain_pdb = os.path.join(output_folder, f'protein_{filename}_{chain_letter}.pdb')
        chain_traj.save(chain_pdb)
        
        # Superpose against reference if provided using PyMOL
        if reference_pdb is not None:
            chain_pdb = align_structures_pymol(chain_pdb, reference_pdb)
        
        # Run pyKVFinder workflow on this chain
        results = pyKVFinder.run_workflow(
            chain_pdb, 
            probe_out=probe_out, 
            volume_cutoff=volume_cutoff, 
            include_depth=True, 
            include_hydropathy=True, 
            ignore_backbone=True
        )
        
        # Prepare output files
        results_toml = os.path.join(output_folder, f'results_{filename}_{chain_letter}.toml')
        output_pdb = os.path.join(output_folder, f'output_{filename}_{chain_letter}.pdb')
        histograms_pdf = os.path.join(output_folder, f'histograms_{filename}_{chain_letter}.pdf')
        cavities_pdb = os.path.join(output_folder, f'cavities_{filename}_{chain_letter}.pdb')
        
        # Export results
        results.export_all(
            fn=results_toml, 
            output=output_pdb, 
            include_frequencies_pdf=True, 
            pdf=histograms_pdf
        )
        
        # Extract cavities information
        areas = results.area
        volumes = results.volume
        max_depth = results.max_depth
        avg_depth = results.avg_depth
        avg_hydropathy = results.avg_hydropathy
        
        # Prepare pocket results DataFrame
        for cavity, area in areas.items():
            all_results.append({
                'Chain': chain_letter,
                'Cavity': cavity, 
                'Volume': volumes.get(cavity),
                'Area': area, 
                'Max_depth': max_depth.get(cavity),
                'Avg_depth': avg_depth.get(cavity), 
                'Avg_hydropathy': avg_hydropathy.get(cavity)
            })
        
        # Write cavities PDB
        pyKVFinder.write_results(
            results_toml, 
            input=chain_pdb, 
            volume=volumes, 
            area=areas, 
            ligand=None, 
            output=cavities_pdb
        )
        
        print(f"\nChain {chain_letter}:")
        print(f"- Protein PDB: {chain_pdb}")
        print(f"- TOML: {results_toml}")
        print(f"- PDB Output: {output_pdb}")
        print(f"- Histograms: {histograms_pdf}")
        print(f"- Cavities PDB: {cavities_pdb}")
    
    # Save combined results to Excel
    excel_file = os.path.join(output_folder, f'pocket_results_{filename}.xlsx')
    pockets_df = pd.DataFrame(all_results)
    pockets_df.to_excel(excel_file, index=False)
    print(f"\nCombined results saved: {excel_file}")

def main():
    parser = argparse.ArgumentParser(description="Process PDB file with pyKVFinder, split by chains")
    parser.add_argument("pdb_file", help="Path to the input PDB file")
    parser.add_argument("--reference", help="Path to reference PDB file for superposition")
    parser.add_argument("--probe_out", type=float, default=8.0, help="Probe out size (default: 8.0)")
    parser.add_argument("--volume_cutoff", type=float, default=50.0, help="Volume cutoff (default: 50.0)")
    args = parser.parse_args()
    
    # Validate files
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file '{args.pdb_file}' does not exist")
        return
    
    if not args.pdb_file.lower().endswith('.pdb'):
        print("Error: Input file must have .pdb extension")
        return
    
    if args.reference and not os.path.exists(args.reference):
        print(f"Error: Reference PDB file '{args.reference}' does not exist")
        return
    
    try:
        process_pdb_file(args.pdb_file, args.probe_out, args.volume_cutoff, args.reference)
    except Exception as e:
        print(f"Error processing PDB file: {str(e)}")

if __name__ == "__main__":
    main()
