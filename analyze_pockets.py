import os
import argparse
import pyKVFinder
import pandas as pd

def process_pdb_file(filepath, probe_out, volume_cutoff):
    """Process PDB file with pyKVFinder."""
    # Ensure output directory exists
    output_folder = 'results'
    os.makedirs(output_folder, exist_ok=True)
    
    # Run pyKVFinder workflow
    results = pyKVFinder.run_workflow(
        filepath, 
        probe_out=probe_out, 
        volume_cutoff=volume_cutoff, 
        include_depth=True, 
        include_hydropathy=True, 
        ignore_backbone=True
    )
    
    # Prepare output files
    filename = os.path.splitext(os.path.basename(filepath))[0]
    results_toml = os.path.join(output_folder, f'results_{filename}.toml')
    output_pdb = os.path.join(output_folder, f'output_{filename}.pdb')
    histograms_pdf = os.path.join(output_folder, f'histograms_{filename}.pdf')
    excel_file = os.path.join(output_folder, f'pocket_results_{filename}.xlsx')
    cavities_pdb = os.path.join(output_folder, f'cavities_{filename}.pdb')
    
    # Export results
    results.export_all(
        fn=results_toml, 
        output=output_pdb, 
        include_frequencies_pdf=True, 
        pdf=histograms_pdf
    )
    
    # Extract cavities information
    cavities = results.cavities
    areas = results.area
    volumes = results.volume
    max_depth = results.max_depth
    avg_depth = results.avg_depth
    avg_hydropathy = results.avg_hydropathy
    
    # Prepare pocket results DataFrame
    pocket_results = []
    for cavity, area in areas.items():
        pocket_results.append({
            'Cavity': cavity, 
            'Volume': volumes.get(cavity),
            'Area': area, 
            'Max_depth': max_depth.get(cavity),
            'Avg_depth': avg_depth.get(cavity), 
            'Avg_hydropathy': avg_hydropathy.get(cavity)
        })
    
    # Convert to DataFrame and save
    pockets_df = pd.DataFrame(pocket_results)
    pockets_df.to_excel(excel_file, index=False)
    
    # Write cavities PDB
    pyKVFinder.write_results(
        results_toml, 
        input=filepath, 
        volume=volumes, 
        area=areas, 
        ligand=None, 
        output=cavities_pdb
    )
    
    print(f"Results saved in {output_folder}:")
    print(f"- TOML: {results_toml}")
    print(f"- PDB Output: {output_pdb}")
    print(f"- Histograms: {histograms_pdf}")
    print(f"- Excel: {excel_file}")
    print(f"- Cavities PDB: {cavities_pdb}")

def main():
    parser = argparse.ArgumentParser(description="Process PDB file with pyKVFinder")
    parser.add_argument("pdb_file", help="Path to the input PDB file")
    parser.add_argument("--probe_out", type=float, default=8.0, help="Probe out size (default: 8.0)")
    parser.add_argument("--volume_cutoff", type=float, default=50.0, help="Volume cutoff (default: 50.0)")
    
    args = parser.parse_args()
    
    # Validate file
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file '{args.pdb_file}' does not exist")
        return
    if not args.pdb_file.lower().endswith('.pdb'):
        print("Error: Input file must have .pdb extension")
        return
    
    try:
        process_pdb_file(args.pdb_file, args.probe_out, args.volume_cutoff)
    except Exception as e:
        print(f"Error processing PDB file: {str(e)}")

if __name__ == "__main__":
    main()
