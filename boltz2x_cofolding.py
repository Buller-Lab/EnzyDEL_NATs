import os
import subprocess
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Process .yml files with boltz predict")
parser.add_argument('--input_folder', required=True, help='Path to input folder containing .yml files')
parser.add_argument('--output_folder', required=True, help='Path to output folder for results')
args = parser.parse_args()

# Get directory paths from arguments
input_directory = args.input_folder
output_base_directory = args.output_folder

# Create output base directory if it doesn't exist
os.makedirs(output_base_directory, exist_ok=True)

# Iterate through all files in the input directory
for filename in os.listdir(input_directory):
    # Check if the file ends with .yml
    if filename.endswith(".yml"):
        # Get the base filename without .yml extension
        base_name = filename.replace(".yml", "")
        
        # Construct the results folder path
        results_folder = os.path.join(output_base_directory, f"botz_results_{base_name}")
        
        # Check if the results folder already exists
        if os.path.exists(results_folder):
            print(f"Skipping {filename}: Results folder '{results_folder}' already exists")
            continue
        
        # Construct the input file path
        input_file = os.path.join(input_directory, filename)
        
        # Construct the command
        command = f"boltz predict {input_file} --use_msa_server --use_potentials --diffusion_samples 30"
        
        # Print the command (optional, for verification)
        print(f"Executing: {command}")
        
        # Execute the command
        try:
            subprocess.run(command, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing command for {filename}: {e}")

print("Processing complete!")
