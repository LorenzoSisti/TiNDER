"""
Script to convert mmCIF files to PDB format.
usage: python cif2pdb.py                   → converte un singolo file
       python cif2pdb.py --batch directory → converte automaticamente tutti i .cif trovati nelle sottocartelle
"""

import sys
import argparse
import logging
import os
from pathlib import Path

from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO


def int_to_chain(i, base=62):
    """(NON MODIFICATO)"""
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
    pass


def rename_chains(structure):
    """(NON MODIFICATO)"""
    next_chain = 0
    chainmap = {c.id: c.id for c in structure.get_chains() if len(c.id) == 1}
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


# 👇 AGGIUNTA: funzione che gestisce la conversione CIF → PDB
def convert_cif_to_pdb(ciffile, pdbfile, verbose=False):
    """Converte un singolo file mmCIF in un file PDB."""
    parser = MMCIFParser()
    # Necessario per Biopython, il nome struttura deve essere max 4 caratteri
    strucid = Path(ciffile).stem.split("_")[1]   # prende "1fdl"

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

    io = PDBIO()
    io.set_structure(structure)
    io.save(str(pdbfile))     # <-- conversione necessaria per Biopython
    print(f"✔ File salvato: {pdbfile}")



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Convert mmCIF to PDB format")
    parser.add_argument("ciffile", nargs="?", help="mmCIF input file")
    parser.add_argument("pdbfile", nargs="?", help="PDB output file. Default based on CIF filename")
    parser.add_argument("--batch", help="Main folder containing the subfolders with CIF files")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode")

    args = parser.parse_args()
    logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG if args.verbose else logging.WARN)

    # ✅ MODALITÀ AUTOMATICA: --batch
    if args.batch:
        main_dir = Path(args.batch)

        if not main_dir.exists():
            sys.exit(f"❌ Directory non trovata: {main_dir}")

        print(f"🔍 Scansione cartelle in: {main_dir}")

        # Scorre tutte le sottocartelle
        for subdir in main_dir.iterdir():
            if subdir.is_dir():
                pdb_id = subdir.name.split("_")[0]  # es: "1fdl_seed1" → "1fdl"

                # cerca file fold_XXXXX_model_X.cif
                for ciffile in subdir.glob("fold_*_model_*.cif"):

                    # Estrai numero model_X dall'ultimo carattere prima del .cif
                    model_num = ciffile.stem.split("_")[-1]  # es: model_0

                    pdb_output = subdir / f"{pdb_id}_{model_num.replace('model_', '')}.pdb"

                    convert_cif_to_pdb(ciffile, pdb_output, verbose=args.verbose)

        sys.exit(0)

    # ✅ Modalità singola (script originale NON modificato)
    if not args.ciffile:
        parser.print_help()
        sys.exit(1)

    ciffile = args.ciffile
    pdbfile = args.pdbfile or ciffile + ".pdb"

    convert_cif_to_pdb(ciffile, pdbfile, verbose=args.verbose)
