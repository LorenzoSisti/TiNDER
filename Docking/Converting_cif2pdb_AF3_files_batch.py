##########################################################################################################
# THIS SCRIPT CONVERTS DOCKED STRUCTURES FROM ALPHAFOLD3 (PROVIDED AS .CIF FILES) INTO STANDARD .PDB FILES.
# IT SUPPORTS BOTH SINGLE-FILE CONVERSION AND AUTOMATIC BATCH PROCESSING ACROSS MULTIPLE SUBDIRECTORIES.

# USAGE EXAMPLES:
#   ▶ python cif2pdb.py input.cif                  → Converts a single CIF file to PDB
#   ▶ python cif2pdb.py --batch /path/to/folders   → Converts all CIF files found in subdirectories

# THIS SCRIPT WAS ADAPTED FROM THE SPENCER BLIVEN'S https://gist.github.com/sbliven/b7cc2c5305aa3652a75a580ae7c6ce33
##########################################################################################################


import sys
import argparse
import logging
import os
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO


def int_to_chain(i, base=62):
    """Convert an integer into a valid 1-character chain ID (A–Z, 0–9, a–z)."""
    if i < 0:
        raise ValueError("positive integers only")
    if base < 0 or 62 < base:
        raise ValueError("Invalid base")

    quot = int(i)//base
    rem = i%base
    if rem < 26:
        letter = chr(ord("A") + rem)
    elif rem < 36:
        letter = str(rem-26)
    else:
        letter = chr(ord("a") + rem - 36)
    if quot == 0:
        return letter
    else:
        return int_to_chain(quot-1, base) + letter


class OutOfChainsError(Exception):
"""Raised when more than 62 unique chain IDs are required (limit exceeded)."""
    pass


def rename_chains(structure):
    """Ensure all chains have valid single-character IDs."""
    next_chain = 0
    chainmap = {c.id: c.id for c in structure.get_chains() if len(c.id) == 1}

# Iterate through all chains and assign new valid IDs where needed
    for o in structure.get_chains():
        if len(o.id) != 1:
            if o.id[0] not in chainmap:
                chainmap[o.id[0]] = o.id
                o.id = o.id[0]
            else:
                c = int_to_chain(next_chain)
                while c in chainmap:
                    next_chain += 1
                    c = int_to_chain(next_chain)
                    if next_chain >= 62:
                        raise OutOfChainsError()
                chainmap[c] = o.id
                o.id = c
    return chainmap


############################################################
# FUNCTION THAT CONVERTS A SINGLE mmCIF FILE TO PDB FORMAT
############################################################
def convert_cif_to_pdb(ciffile, pdbfile, verbose=False):
    """Converte un singolo file mmCIF in un file PDB."""
    parser = MMCIFParser()
    # Biopython requires the structure ID to be short (≤4 characters)
    strucid = Path(ciffile).stem.split("_")[1]   # Extract ID such as "1fdl"

    structure = parser.get_structure(strucid, ciffile)

    try:
        chainmap = rename_chains(structure)
    except OutOfChainsError:
        logging.error(f"❌ Troppi chain ID nel file {ciffile}")
        return

    if verbose:
        logging.info(f"✅ Converting: {ciffile}")
        for new, old in chainmap.items():
            if new != old:
                logging.info(f"   Chain renamed {old} → {new}")

    # Save structure as PDB
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdbfile))     # <-- Biopython requires a string path
    print(f"✔ File salvato: {pdbfile}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert mmCIF to PDB format")
    parser.add_argument("ciffile", nargs="?", help="mmCIF input file")
    parser.add_argument("pdbfile", nargs="?", help="PDB output file. Default based on CIF filename")
    parser.add_argument("--batch", help="Main folder containing the subfolders with CIF files")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")

    args = parser.parse_args()
    logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG if args.verbose else logging.WARN)

    ############################################################
    # AUTOMATIC MODE (--batch): CONVERT ALL .CIF FILES IN SUBDIRECTORIES
    ############################################################
    if args.batch:
        main_dir = Path(args.batch)

        if not main_dir.exists():
            sys.exit(f"❌ Directory non trovata: {main_dir}")

        print(f"🔍 Scansione cartelle in: {main_dir}")

        # Loop through all subdirectories
        for subdir in main_dir.iterdir():
            if subdir.is_dir():
                pdb_id = subdir.name.split("_")[0]  # es: "1fdl_seed1" → "1fdl"

                # Search for CIF files like "fold_XXXXX_model_X.cif" (AF3 output)
                for ciffile in subdir.glob("fold_*_model_*.cif"):

                    # Extract model number (e.g., model_0)
                    model_num = ciffile.stem.split("_")[-1]  # es: model_0
                    # Create output filename: e.g., 1fdl_0.pdb
                    pdb_output = subdir / f"{pdb_id}_{model_num.replace('model_', '')}.pdb"

                    convert_cif_to_pdb(ciffile, pdb_output, verbose=args.verbose)

        sys.exit(0)

    ############################################################
    # SINGLE-FILE MODE: CONVERT ONE CIF → PDB
    ############################################################
    if not args.ciffile:
        parser.print_help()
        sys.exit(1)

    ciffile = args.ciffile
    pdbfile = args.pdbfile or ciffile + ".pdb"

    convert_cif_to_pdb(ciffile, pdbfile, verbose=args.verbose)
