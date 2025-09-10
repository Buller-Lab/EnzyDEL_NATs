import os
import subprocess
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Process .yml files with boltz predict")
parser.add_argument('--input_folder', required=True, help='Path to input folder containing .yml files')
args = parser.parse_args()

# Get directory path from arguments
input_directory = args.input_folder

# Iterate through all files in the input directory
for filename in os.listdir(input_directory):
    # Check if the file ends with .yml
    if filename.endswith(".yml"):
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
