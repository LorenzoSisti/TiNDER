# THIS SCRIPT WAS USED IN ORDER TO:

# CALCULATE THE RSA FOR EACH RESIDUE OF THE AB::AG COMPLEXES TO BE DOCKED, IN ORDER TO ELIMINATE THE PREDICTED PARATOPE RESIDUES (proABC-2) TO BE REMOVED
# CALCULATE THE RSA FOR ALL THE OVER 2000 AB::AG COMPLEXES, IN ORDER TO INVESTIGATE INTERFACE PROPENSITY, NON-INTERACTING SURFACES PROPENSITY AND CORE PROPENSITY

# Importing libraries
import freesasa
import pandas as pd
import math
from pathlib import Path

# Defining the directories with the PDB file collection, and where to save the output .csv files
pdb_folder = Path('/Users/lorenzosisti/Downloads/database_settembre_renamed/')
output_folder = Path('/Users/lorenzosisti/Downloads/SASA_for_ab_ag_int_characterization/')

# Creating a dictionary containing the 20 experimentally derived maximum allowed solvent accessibility values reported in Tien et al. (2023).
max_acc = {
    'ALA': 121.0,
    'ARG': 265.0,
    'ASN': 187.0,
    'ASP': 187.0,
    'CYS': 148.0,
    'GLN': 214.0,
    'GLU': 214.0,
    'GLY': 97.0,
    'HIS': 216.0,
    'ILE': 195.0,
    'LEU': 191.0,
    'LYS': 230.0,
    'MET': 203.0,
    'PHE': 228.0,
    'PRO': 154.0,
    'SER': 143.0,
    'THR': 163.0,
    'TRP': 264.0,
    'TYR': 255.0,
    'VAL': 165.0
}

# Defining a function to compute the relative solvent accessibility (RSA) of each residue
def compute_rsa(row, max_acc_dict):

    aa = row["ResName"]
    sasa = row["SASA"]

    if aa in max_acc_dict and sasa is not None:
        maxval = max_acc_dict[aa]
        if maxval and not math.isclose(maxval, 0.0):
            return sasa / maxval

    return float("nan")

# Function to process each PDB file (i.e. computing SASA and RSA)
def process_pdb(pdb_path, max_acc_dict):

    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)

    records = []
    for i in range(structure.nAtoms()):
        chain = structure.chainLabel(i)
        resnum = structure.residueNumber(i)
        resname = structure.residueName(i).upper()
        sasa_atom = result.atomArea(i)
        records.append([chain, resnum, resname, sasa_atom])

    df_atoms = pd.DataFrame(records, columns=["Chain", "ResNum", "ResName", "AtomSASA"])
    df_res = df_atoms.groupby(["Chain", "ResNum", "ResName"], as_index=False)["AtomSASA"].sum()
    df_res = df_res.rename(columns={"AtomSASA": "SASA"})

    df_res["RSA"] = df_res.apply(lambda row: compute_rsa(row, max_acc_dict), axis=1)

    df_res["PDB"] = Path(pdb_path).name

    return df_res

# Batch run on all the PDB in the folder of interest
pdb_files = sorted([p for p in pdb_folder.iterdir() if p.suffix.lower() in (".pdb", ".ent")])
if not pdb_files:
    raise FileNotFoundError(f"Nessun file .pdb trovato in {pdb_folder}")

combined_list = []
failed = []
combined_csv = output_folder / "combined_sasa_rsa.csv"

for pdb in pdb_files:
    try:
        print(f"Processing {pdb.name} ...")
        df_res = process_pdb(pdb, max_acc)
        out_csv = output_folder / f"{pdb.stem}_sasa_rsa.csv"
        df_res.to_csv(out_csv, index=False)
        print(f"  → salvato: {out_csv}")
        combined_list.append(df_res)
    except Exception as e:
        print(f"  ! Errore su {pdb.name}: {e}")
        failed.append((pdb.name, str(e)))

if combined_list:
    combined_df = pd.concat(combined_list, ignore_index=True)
    combined_df.to_csv(combined_csv, index=False)
    print(f"\n✅ File combinato salvato: {combined_csv} (righe: {len(combined_df)})")
else:
    print("\n⚠️ Nessun file processato correttamente, nessun file combinato creato.")

if failed:
    print("\nAlcuni file non sono stati processati:")
    for name, err in failed:
        print(" -", name, ":", err)
else:
    print("\nTutti i file processati con successo.")