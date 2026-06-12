import os

base_dir = "/Users/lorenzosisti/Downloads/AF3_docking_50_complexes_11_06_2026"

for folder_name in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder_name)

    if not os.path.isdir(folder_path):
        continue
    if not folder_name.endswith("_seed2"):
        continue

    print(f"\nCartella: {folder_name}")

    for filename in sorted(os.listdir(folder_path)):
        if not filename.endswith("_renamed.pdb"):
            continue

        # filename è tipo "2r6p_0_renamed.pdb"
        parts = filename.split("_")   # ["2r6p", "0", "renamed.pdb"]
        number = int(parts[1])
        parts[1] = str(number + 5)
        new_filename = "_".join(parts)

        os.rename(
            os.path.join(folder_path, filename),
            os.path.join(folder_path, new_filename)
        )
        print(f"  {filename}  ->  {new_filename}")

print("\nDone.")