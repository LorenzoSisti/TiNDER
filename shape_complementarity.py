import subprocess
import os

# Set the paths for the Rosetta binary and XML script
rosetta_binary = "/path/to/rosetta/source/bin/rosetta_scripts.default.macosclangrelease"
xml_file = "/path/to/your/xml_script/shape_complementarity.xml"

input_folder = "/path/to/your/non/redundant/pdb/files/directory/" # Here, we have all the PDBs from the database filtering process
output_score_folder = "/path/to/the/directory/where/you/want/to/save/the/scores/"

log_file_path = os.path.join(output_score_folder, "output_rosetta_shape_log.txt")

# Clean the log file (if it already exists)
with open(log_file_path, "w") as log_file:
    log_file.write("=== LOG ESECUZIONE ROSETTA ===\n\n")

# Get all the PDBs from the input folder
pdb_files = [f for f in os.listdir(input_folder) if f.endswith(".pdb")]

if not pdb_files:
    print("No PDB files found in the folder.")
    exit()

# Iterate on each PDB file
for pdb_file in pdb_files:
    input_pdb_path = os.path.join(input_folder, pdb_file)

    # Remove the .pdb extension
    pdb_id = os.path.splitext(pdb_file)[0]  # es: "1a3r"
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
    result = subprocess.run(command, capture_output=True, text=True)

    # Write the log file
    with open(log_file_path, "a") as log_file:
        log_file.write(f"\n=== Output per {pdb_file} ===\n")
        log_file.write("--- STDOUT ---\n")
        log_file.write(result.stdout)
        log_file.write("\n--- STDERR ---\n")
        log_file.write(result.stderr)
        log_file.write(f"\nReturn code: {result.returncode}\n")
        log_file.write("="*40 + "\n")

    # Check if the command was successful
    if result.returncode != 0:
        print(f"❌ Error during execution on {pdb_file}.")
    else:
        print(f"✅ Success on {pdb_file}.")


### To merge all the score files in a single file ###

# Get all the score files .sc
score_files = [f for f in os.listdir(output_score_folder) if f.startswith("score_") and f.endswith(".sc")]
score_totale_path = os.path.join(output_score_folder, "score_totale.sc")

with open(score_totale_path, "w") as outfile:
    for i, score_file in enumerate(score_files):
        score_path = os.path.join(output_score_folder, score_file)
        with open(score_path, "r") as infile:
            lines = infile.readlines()

            # If it's the first file, write all lines (including header)
            # If it's not the first file, skip the header
            if i == 0:
                outfile.writelines(lines)
            else:
                outfile.writelines(lines[1:])  # Skip header