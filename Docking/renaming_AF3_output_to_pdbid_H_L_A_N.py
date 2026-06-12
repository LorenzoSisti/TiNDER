import os

folder = "/Users/lorenzosisti/Downloads/docked_structures_renamed_AF3_11_06"

for filename in sorted(os.listdir(folder)):
    if not filename.endswith("_renamed.pdb"):
        continue

    # "1hez_0_renamed.pdb" -> ["1hez", "0", "renamed.pdb"]
    parts = filename.split("_")
    pdbid = parts[0]   # es. "1hez"
    number = parts[1]  # es. "0"

    new_filename = f"{pdbid}_H_L_A_{number}.pdb"

    os.rename(
        os.path.join(folder, filename),
        os.path.join(folder, new_filename)
    )
    print(f"  {filename}  ->  {new_filename}")

print("\nDone.")