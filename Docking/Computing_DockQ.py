##########################################################################################################
# PURPOSE OF THIS SCRIPT:
# THIS SCRIPT AUTOMATES THE CALCULATION OF DOCKQ SCORES FOR ANTIBODY-ANTIGEN COMPLEXES.
# IT MERGES HEAVY (H) AND LIGHT (L) ANTIBODY CHAINS INTO A SINGLE CHAIN (CHAIN A),
# RENAMES THE ANTIGEN AS CHAIN B, AND THEN COMPUTES INTERFACE QUALITY METRICS
# (DOCKQ, IRMSD, LRMSD, FNAT, AND RMSD) BY COMPARING PREDICTED MODELS TO NATIVE STRUCTURES.

# THIS SCRIPT TAKE AS INPUT NATIVE/REFERENCE STRUCTURES FILES NAMED PDBID_H_L_Ag.pdb where: 
    # H is (independently from the letter) the heavy chain
    # L is (independently from the letter) the light chain
    # Ag is the antigen chain
# THIS SCRIPT TAKE AS INPUT DOCKED STRUCTURES FILES NAMED PDBID_H_L_Ag_N.pdb where: 
    # H, L, Ag are as defined before
    # N is a number identifying the docking pose (as more than one pose was generated per PDB entry
##########################################################################################################

import os
import tempfile
import string
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain, Residue
from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
from Bio.PDB.Superimposer import Superimposer

############################################################
# FUNCTION: MERGE ANTIBODY HEAVY (H) + LIGHT (L) CHAINS INTO ONE (CHAIN A) AND THE ANTIGEN INTO CHAIN B.
############################################################

def merge_chains_and_rename_in_memory_safe(input_pdb, chain_H, chain_L, chain_Ag):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("src", input_pdb)
    model = next(structure.get_models())

    # Create a new clean structure and model to avoid residue ID conflicts
    new_structure = Structure.Structure("merged")
    new_model = Model.Model(0)
    new_structure.add(new_model)

    chainA = Chain.Chain("A")
    new_model.add(chainA)

    new_resid = 1 # Start renumbering residues from 1

    # Merge heavy and light chains into chain A
    for cid in (chain_H, chain_L):
        if cid not in model:
            continue

        for old_res in model[cid]:
            # Create a new residue with sequential numbering
            new_res = Residue.Residue((" ", new_resid, " "), old_res.get_resname(), old_res.segid)

            # Copy all atoms from old residue to new one
            for atom in old_res.get_atoms():
                new_atom = atom.copy()
                if new_atom.occupancy is None:
                    new_atom.set_occupancy(1.00)
                new_res.add(new_atom)

            chainA.add(new_res)      
            new_resid += 1           # Increment residue number for next residue

    # Create chain B for antigen
    if chain_Ag in model:
        chainB = Chain.Chain("B")
        new_model.add(chainB)

        for old_res in model[chain_Ag]:
            # Keep original residue numbering for antigen
            new_res = Residue.Residue(old_res.id, old_res.get_resname(), old_res.segid)
            for atom in old_res.get_atoms():
                new_atom = atom.copy()
                if new_atom.occupancy is None:
                    new_atom.set_occupancy(1.00)
                new_res.add(new_atom)
            chainB.add(new_res)

    # Write the clean structure to a temporary file for DockQ analysis
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
    tmp_path = tmp.name
    tmp.close()

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(tmp_path)

    return tmp_path



############################################################
# RMSD (for debugging purposes)
############################################################

def rmsd_biopython(model_path, native_path):
    
    """Compute RMSD between model and native structures using Cα atoms."""
    
    parser = PDBParser(QUIET=True)
    model = parser.get_structure("model", model_path)
    native = parser.get_structure("native", native_path)

    # Extract only alpha-carbon atoms for RMSD calculation
    model_atoms = [a for a in model.get_atoms() if a.get_name() == "CA"]
    native_atoms = [a for a in native.get_atoms() if a.get_name() == "CA"]

    n = min(len(model_atoms), len(native_atoms))
    if n == 0:
        return float("nan")

    # Superimpose model onto native structure
    sup = Superimposer()
    sup.set_atoms(native_atoms[:n], model_atoms[:n])
    sup.apply(model.get_atoms())

    return sup.rms


############################################################
# MAIN LOOP — RUN DOCKQ COMPARISONS FOR EACH MODEL
############################################################

native_dir = "/path/to/native/structures"
model_dir = "/path/to/docked/structures"

out_file = "DockQ_results.csv"

# Create output CSV 
with open(out_file, "w") as f:
    f.write("PDB_ID,Model,ChainHeavy,ChainLight,ChainAntigen,DockQ,iRMSD,LRMSD,fnat,RMSD\n")

# Loop over all model PDBs in the directory
for model_file in sorted(os.listdir(model_dir)):
    if not model_file.endswith(".pdb"):
        continue

    print(f"\n➡️ Calculating DockQ for model: {model_file}")

    # Parse naming convention: PDBID_H_L_Ag.pdb
    parts = model_file.replace(".pdb", "").split("_")
    pdb_id, chain_H, chain_L, chain_Ag = parts[:4]

    model_path = os.path.join(model_dir, model_file)
    native_path = os.path.join(native_dir, f"{pdb_id}_{chain_H}_{chain_L}_{chain_Ag}.pdb")
    
    # Skip if native reference is missing
    if not os.path.exists(native_path):
        print(f"❌ File nativo mancante: {native_path}")
        continue

    model_tmp = native_tmp = None
    try:
        
        # Generate clean temporary merged PDBs for both model and native
        model_tmp = merge_chains_and_rename_in_memory_safe(model_path, chain_H, chain_L, chain_Ag)
        native_tmp = merge_chains_and_rename_in_memory_safe(native_path, chain_H, chain_L, chain_Ag)

        model = load_PDB(model_tmp)
        native = load_PDB(native_tmp)
        
        # Run DockQ on the antibody-antigen interface (A–B)
        results = run_on_all_native_interfaces(model, native, chain_map={"A": "A", "B": "B"})
        interface = list(results[0].keys())[0]
        
        # Extract key DockQ metrics
        dockq = results[0][interface]["DockQ"]
        irmsd = results[0][interface]["iRMSD"]
        lrmsd = results[0][interface]["LRMSD"]
        fnat = results[0][interface]["fnat"]
        
        # Compute RMSD for debugging
        rmsd_debug = rmsd_biopython(model_tmp, native_tmp)
        
        # Append results to the CSV file
        with open(out_file, "a") as f:
            f.write(f"{pdb_id},{model_file},{chain_H},{chain_L},{chain_Ag},"
                    f"{dockq},{irmsd},{lrmsd},{fnat},{rmsd_debug}\n")

        print(f"✅ DockQ = {dockq:.3f}")

    except Exception as e:
        print(f"❌ Errore su {model_file}: {e}")

    finally:        
        
        # Clean up temporary files
        if model_tmp and os.path.exists(model_tmp): os.remove(model_tmp)
        if native_tmp and os.path.exists(native_tmp): os.remove(native_tmp)

