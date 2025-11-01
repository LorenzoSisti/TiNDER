# Docking and docking pose evaluation

Descrizione generale

## Structure pre-processing

Download PDB file from Protein Data Bank, poi aggiustali con [PDBfixer](https://github.com/openmm/pdbfixer) tramite:

```
pdbfixer PDB_ID.pdb --replace-nonstandard --add-residues --keep-heterogens=all --output=PDB_ID_fixed.pdb
```

Work in progress...
