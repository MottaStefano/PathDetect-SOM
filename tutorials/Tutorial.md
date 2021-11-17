# PathDetect-SOM Tutorial: 

In this tutorial, you will learn how to analyize ligand binding molecular dynamics (MD) simulations with the PathDetect-SOM tool.
In the SIMS folder you will find 15 replicas of steered MD simulations of a ligand unbinding from the Hypoxia inducible Factor-2α (HIF-2α) protein[^1] that are a small subset of the original data. 

1. ## **Generate coordinate file**

As a first step, coordinates for the atoms of interest should be extracted from the simulation. This require the construction of an index file using the gmx make_ndx tool, creating a group that contains the atoms of interest. Selected atoms should describe both the binding site and the mouth at the entrance of the binding site. Ideally, both atoms from backbone and from large or polar/charged sidechains should be included when the side chain dynamics and interactions are relevant for binding. Similarly, selected ligand atoms should well describe the core molecular structure and all the relevant lateral groups. 
The index file for the case here presented is already present in the SIMS folder. Then the atoms coordinates can be extracted using the gromacs command gmx traj:

    mkdir COORDS
    for REP in {01..15}; do
    gmx traj -f SIMS/Rep_${REP}.xtc -s SIMS/Protein-lig.gro -ox COORDS/coord-${REP}.xvg -n SIMS/index.ndx <<EOC
    16
    EOC
    done

Here we used a bash loop to perform the calculation over all the 15 replicas. The output xvg files that contains the coordinates for the atoms of interest are saved in the COORDS folder.

2. ## **Train the SOM**

In the second step, we will train a Self-Organizing Map using all the intermolecular distances between selected protein and ligand atoms. This can be done using the PathDetect-SOM script. To obtain an help page from the tool simply run:

    PathDetect-SOM.r

or

    PathDetect-SOM.r --help
    




[^1]: Callea, L.; Bonati, L.; Motta, S. Metadynamics-Based Approaches for Modeling the Hypoxia-Inducible Factor 2α Ligand Binding Process. Journal of Chemical Theory and Computation 2021, 18. https://doi.org/10.1021/acs.jctc.1c00114.
