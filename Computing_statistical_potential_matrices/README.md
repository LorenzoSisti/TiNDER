# Statistical potentials

This folder contains the script necessary to compute whole-interface and radially distributed symmetric and asymmetric pairwise potentials (or PMF, potential of mean force) from the antibody-antigen complexes representative collection, and the scripts to assign to each predicted complex the corresponding statistical potential value.

## Computing pairwise potentials

The `Computing_whole_interface_pairwise_potentials.R` script computes coarse grained whole-interface symmetric and asymmetric PMF from a set of antibody-antigen complexes representative PDB files collection. 

The `Computing_radially_distributed_pairwise_potentials.R` script computes coarse grained radially distributed symmetric and asymmetric PMF from a set of antibody-antigen complexes representative PDB files collection. Its output is a wide-format .csv file ...

## Assigning statistical potentials to predicted complexes

The `Assigning_whole_interface_statistical_potentials_to_docked_structures.R` script uses as input the .csv files obtained via the `Computing_whole_interface_pairwise_potentials.R` script, and sums the whole-interface pairwise interaction potentials to assign to each predicted complex a statistical potential value. 


The `Assigning_g_r_potentials_for_docking_pose_opt.R.R` script uses as input the wide .csv files obtained via the `Computing_radially_distributed_pairwise_potentials.R` script, and sums the stratified pairwise interaction potentials to assign to each predicted complex a more *fine-grained* statistical potential value. 

