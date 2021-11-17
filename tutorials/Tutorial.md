# PathDetect-SOM Tutorial: 

In this tutorial, you will learn how to analyize ligand binding molecular dynamics (MD) simulations with the PathDetect-SOM tool.
In the SIMS folder you will find 15 replicas of steered MD simulations of a ligand unbinding from the Hypoxia inducible Factor-2α (HIF-2α) protein[^1] that are a small subset of the original data. 

1. ## **Generate coordinate file**

As a first step, coordinates for the atoms of interest should be extracted from the simulation. This require the construction of an index file using the gmx make_ndx tool, creating a group that contains the atoms of interest. The index file for the case here presented is already present in the SIMS folder. Then the atoms coordinates can be extracted
using the gromacs command gmx traj:

    mkdir COORDS
    for REP in {01..15}; do
    gmx traj -f SIMS/Rep_${REP}.xtc -s SIMS/Protein-lig.gro -ox COORDS/coord-${REP}.xvg -n SIMS/index.ndx <<EOC
    16
    EOC
    done

Here we used a bash loop to perform the calculation over all the 15 replicas. The output xvg files that contains the coordinates for the atoms of interest are saved in the COORDS folder.







[^1]: Callea, L.; Bonati, L.; Motta, S. Metadynamics-Based Approaches for Modeling the Hypoxia-Inducible Factor 2α Ligand Binding Process. Journal of Chemical Theory and Computation 2021, 18. https://doi.org/10.1021/acs.jctc.1c00114.
