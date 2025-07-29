### Script that is used to run the .xml ROSETTASCRIPT to calculate the electrostatic complementarity at the Ab:::Ag interface ###
### https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/xsd/simple_metric_ElectrostaticComplementarityMetric_type ###

import subprocess
import os

# Set the paths for the Rosetta binary and XML script
rosetta_binary = "/path/to/rosetta/source/bin/rosetta_scripts.default.macosclangrelease"
xml_file = "/path/to/your/xml_script/electrostatic_complementarity.xml"

input_folder = "/path/to/your/non/redundant/pdb/files/directory/" # Here, we have all the PDBs from the database filtering process
output_score_folder = "/path/to/the/directory/where/you/want/to/save/the/scores/"

log_file_path = os.path.join(output_score_folder, "output_rosetta_electrostatic_log.txt")

# Clean the log file (if it already exists)
with open(log_file_path, "w") as log_file:
    log_file.write("=== LOG ESECUZIONE ROSETTA (ELECTROSTATIC) ===\n\n")

# Get all the PDBs from the input folder
pdb_files = [f for f in os.listdir(input_folder) if f.endswith(".pdb")]

if not pdb_files:
    print("❌ No PDB files found in the folder.")
    exit()

# Iterate on each PDB file
for pdb_file in pdb_files:
    input_pdb_path = os.path.join(input_folder, pdb_file)

    # Remove the .pdb extension
    pdb_id = os.path.splitext(pdb_file)[0]  # es. 1a3r
    score_filename = f"score_{pdb_id}.sc"
    score_path = os.path.join(output_score_folder, score_filename)

    command = [
        rosetta_binary,
        "-parser:protocol", xml_file,
        "-in:file:s", input_pdb_path,
        "-out:file:scorefile", score_path,
        "-out:no_nstruct_label",
        "-nstruct", "1",
        "-overwrite"
    ]

    print(f"Executing Rosetta on: {pdb_file}")
    result = subprocess.run(command, capture_output=True, text=True, cwd=output_score_folder)


    # Write the log file
    with open(log_file_path, "a") as log_file:
        log_file.write(f"\n=== Output for {pdb_file} ===\n")
        log_file.write("--- STDOUT ---\n")
        log_file.write(result.stdout)
        log_file.write("\n--- STDERR ---\n")
        log_file.write(result.stderr)
        log_file.write(f"\nReturn code: {result.returncode}\n")
        log_file.write("="*40 + "\n")

    if result.returncode != 0:
        print(f"❌ Error on {pdb_file}")
    else:
        print(f"✅ Success on {pdb_file}")

### To merge all the score files in a single file ###

score_files = [f for f in os.listdir(output_score_folder) if f.startswith("score_") and f.endswith(".sc")]
score_totale_path = os.path.join(output_score_folder, "score_totale_electrostatic.sc")

with open(score_totale_path, "w") as outfile:
    for i, score_file in enumerate(score_files):
        score_path = os.path.join(output_score_folder, score_file)
        with open(score_path, "r") as infile:
            lines = infile.readlines()
            if i == 0:
                outfile.writelines(lines)
            else:
                outfile.writelines(lines[1:])  # skip header

# Extract only the description (PDB name) and ec_avg column and save to CSV
import pandas as pd

score_totale_path = os.path.join(output_score_folder, "score_totale_electrostatic.sc")
output_filtered_path = os.path.join(output_score_folder, "ec_avg_per_structure.csv")

# Read the score file
with open(score_totale_path, "r") as f:
    lines = [line.strip() for line in f if line.startswith("SCORE:")]

# Use the second line as header and the rest as data
header_line = lines[1].split()
data_lines = [line.split() for line in lines[2:]]  

# Create a DataFrame
df = pd.DataFrame(data_lines, columns=header_line)

# Save ec_avg as float 
df["ec_avg"] = pd.to_numeric(df["ec_avg"], errors="coerce")

# Extract only the description and ec_avg columns
df_filtered = df[["description", "ec_avg"]]

# Save to CSV
df_filtered.to_csv(output_filtered_path, index=False)

print(f"\n✅ File saved as: {output_filtered_path}")
