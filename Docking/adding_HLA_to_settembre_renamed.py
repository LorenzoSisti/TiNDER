import os

folder = "/Users/lorenzosisti/Downloads/database_settembre_renamed_copy"

for filename in sorted(os.listdir(folder)):
    if not filename.endswith(".pdb"):
        continue

    pdbid = filename.replace(".pdb", "")   # "1hez.pdb" -> "1hez"
    new_filename = f"{pdbid}_H_L_A.pdb"

    os.rename(
        os.path.join(folder, filename),
        os.path.join(folder, new_filename)
    )
    print(f"  {filename}  ->  {new_filename}")

print("\nDone.")