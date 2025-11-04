import os
import tempfile
import string
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
from Bio.PDB.Superimposer import Superimposer

############################################################
# FUNZIONE: merge H+L in A (rinumerati), antigene in B
# ✅ senza warning, perché creiamo una NUOVA struttura pulita
############################################################

def merge_chains_and_rename_in_memory_safe(input_pdb, chain_H, chain_L, chain_Ag):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("src", input_pdb)
    model = next(structure.get_models())

    new_structure = Structure.Structure("merged")
    new_model = Model.Model(0)
    new_structure.add(new_model)

    chainA = Chain.Chain("A")
    new_model.add(chainA)

    new_resid = 1

    # Merge catene H + L → catena A rinumerata
    for cid in (chain_H, chain_L):
        if cid not in model:
            continue

        for old_res in model[cid]:
            new_res = Residue.Residue((" ", new_resid, " "), old_res.get_resname(), old_res.segid)

            for atom in old_res.get_atoms():
                new_atom = atom.copy()
                if new_atom.occupancy is None:
                    new_atom.set_occupancy(1.00)
                new_res.add(new_atom)

            chainA.add(new_res)      # ✅ dentro il ciclo
            new_resid += 1           # ✅ dentro il ciclo

    # Chain B = antigene
    if chain_Ag in model:
        chainB = Chain.Chain("B")
        new_model.add(chainB)

        for old_res in model[chain_Ag]:
            new_res = Residue.Residue(old_res.id, old_res.get_resname(), old_res.segid)
            for atom in old_res.get_atoms():
                new_atom = atom.copy()
                if new_atom.occupancy is None:
                    new_atom.set_occupancy(1.00)
                new_res.add(new_atom)
            chainB.add(new_res)

    # Scrive struttura su file temporaneo per DockQ
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    tmp_path = tmp.name
    tmp.close()

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(tmp_path)

    return tmp_path



############################################################
# RMSD (solo per debug)
############################################################

def rmsd_biopython(model_path, native_path):
    parser = PDBParser(QUIET=True)
    model = parser.get_structure("model", model_path)
    native = parser.get_structure("native", native_path)

    model_atoms = [a for a in model.get_atoms() if a.get_name() == "CA"]
    native_atoms = [a for a in native.get_atoms() if a.get_name() == "CA"]

    n = min(len(model_atoms), len(native_atoms))
    if n == 0:
        return float("nan")

    sup = Superimposer()
    sup.set_atoms(native_atoms[:n], model_atoms[:n])
    sup.apply(model.get_atoms())

    return sup.rms


############################################################
# MAIN LOOP
############################################################

native_dir = "/Users/lorenzosisti/Downloads/shared"
model_dir = "/Users/lorenzosisti/Downloads/models"

out_file = "DockQ_results.csv"
with open(out_file, "w") as f:
    f.write("PDB_ID,Model,ChainHeavy,ChainLight,ChainAntigen,DockQ,iRMSD,LRMSD,fnat,RMSD\n")

for model_file in sorted(os.listdir(model_dir)):
    if not model_file.endswith(".pdb"):
        continue

    print(f"\n➡️ Calcolo DockQ per modello: {model_file}")

    parts = model_file.replace(".pdb", "").split("_")
    pdb_id, chain_H, chain_L, chain_Ag = parts[:4]

    model_path = os.path.join(model_dir, model_file)
    native_path = os.path.join(native_dir, f"{pdb_id}_{chain_H}_{chain_L}_{chain_Ag}.pdb")

    if not os.path.exists(native_path):
        print(f"❌ File nativo mancante: {native_path}")
        continue

    model_tmp = native_tmp = None
    try:
        model_tmp = merge_chains_and_rename_in_memory_safe(model_path, chain_H, chain_L, chain_Ag)
        native_tmp = merge_chains_and_rename_in_memory_safe(native_path, chain_H, chain_L, chain_Ag)

        model = load_PDB(model_tmp)
        native = load_PDB(native_tmp)

        results = run_on_all_native_interfaces(model, native, chain_map={"A": "A", "B": "B"})
        interface = list(results[0].keys())[0]

        dockq = results[0][interface]["DockQ"]
        irmsd = results[0][interface]["iRMSD"]
        lrmsd = results[0][interface]["LRMSD"]
        fnat = results[0][interface]["fnat"]

        rmsd_debug = rmsd_biopython(model_tmp, native_tmp)

        with open(out_file, "a") as f:
            f.write(f"{pdb_id},{model_file},{chain_H},{chain_L},{chain_Ag},"
                    f"{dockq},{irmsd},{lrmsd},{fnat},{rmsd_debug}\n")

        print(f"✅ DockQ = {dockq:.3f}")

    except Exception as e:
        print(f"❌ Errore su {model_file}: {e}")

    finally:
        if model_tmp and os.path.exists(model_tmp): os.remove(model_tmp)
        if native_tmp and os.path.exists(native_tmp): os.remove(native_tmp)

