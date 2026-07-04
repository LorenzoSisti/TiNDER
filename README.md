# M.Sc. thesis internship
TiNDER (acronym for **T**hesis **i**nternship: **N**ative-like **D**ocking poses **E**valuation & **R**anking) is a protocol for scoring near-native docking poses through *radially distributed* statistical potentials, developed during my MSc thesis internship in Molecular Computational Biophysics at the [Center for Life Nano- & Neuro-Science](https://www.iit.it/it/clns-sapienza) of the Italian Institute of Technology (IIT), under the supervision of [Edoardo Milanetti](https://scholar.google.it/citations?user=Pc8OAWsAAAAJ&hl=it).

## Abstract

*Work in progress*

## Git repository organization

The `dummy_sabdab_database_filtering.R` is a dummy (and simplified) version of the script used to obtain the antibody-antigen complexes dataset used for deriving the scoring function and for the docking simulations. Starting from the entire antibody-antigen complexes collection present in the Protein Data Bank and downloaded from [SAbDab](https://sabdab.opig.stats.ox.ac.uk), renumbered according to the *Chothia* scheme, this script retains exclusively trimeric proteic complexes, that would be later subjected to MSA via the CD-HIT pipeline to keep non-redundant complexes, and to energy minimization via GROMACS to avoid non-physical contacts.

### Antibody_antigen_interfaces_characterization

This repository contains the scripts used to characterise each antibody-antigen interface. The characterisation is necessary not only to quantify the aminoacids heterogeneity at the binding region, but also to extract a smaller subset of complexes representative of the interface diversity to be subjected to docking via the AlphaFold3 Webserver.

### knowledge_based_statistical_potentials

This directory contains:

- The computational pipeline for computing symmetric and asymmetric whole-interface statistical potentials, which is the gold standad currently used in literature.
- The computational pipeline for computing two kinds of symmetric and asymmetric regional statistical potentials, which account for the compositional gradients present at the interface.
- The computational pipeline necessary to assign a statistical potential value to each docking pose.

### docking

### Statistical_potential_performances_saturation
