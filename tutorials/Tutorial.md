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

Here we used a bash loop to perform the calculation over all the 15 replicas. The output xvg files that contains the coordinates for the atoms of interest are saved in the COORDS folder. Note that the PathDetect-SOM tool will read the trajectoris in alphabetical order, so it is recomended to name the files with the same name followed by numbers with leading zeros.

2. ## **Train the SOM**

In the second step, we will train a Self-Organizing Map using all the intermolecular distances between selected protein and ligand atoms. This can be done using the PathDetect-SOM script. To obtain an help page from the tool simply run:

    PathDetect-SOM.r
or
    PathDetect-SOM.r --help
    
The list of the available options will be printed. The mandatory arguments are --folder and --rep, with which one can indicate the folder for the xvg coordinate file and the number of replicas of the simulation. Other options that are usefull for the present case are:
- --type       - In which we specify the type of feature to be used. In this case we want to train the SOM with distances among atoms and not directly with coordinates, so we will choose dRMSD.
- --lig        - that specify that we are treating a ligand-protein system, and thus only intermolecular distances will be computed. Here the range of ligand atoms within the selection should be specified. In this case we have selected a total of 45 atoms, with the fist 40 atoms belonging to protein, and the last five to the ligand. For this reason here we specify the ligand atom range: 41-45.
- --cutoff      - that specify the capping value in nanometers for the distances.
- --path_clus   - that specify the type of pathway clustering. Here we use a time-dependent clustering because in steered MD the replicas evolves in parallel (at the same speed).
- --out         - the name of the output folder (that will be automatically created by the tool).
- --seed        - The random seed for the calculation. This is usefull to ensure reproducibilty.

To run the analysis:

    PathDetect-SOM.r --folder=COORDS --rep=15 --type=dRMSD --lig=41-45 --cutoff=1.2 --path_clus=dependent --out=SOM --seed=30076

Due to the small dataset, the calculation would require just few minutes (depending on the hardware).
At the end of the calculation within the folder COORDS, you will find other sub-folders containing files usefull for the analysis. 

3. ## **Decide the number of cluster**

By default the SOM neurons are further clustered in 10 geometric clusters that represents macro-conformations of the system. However this value may not be the optimal number of cluster. To help users in deciding the optimal number of clusters, results for silhouette analysis are reported within the folder Clusters/Silhouette. For the case of the tutorial a good number of cluster could be 13 (a local maximum in the Silhouette profile (see file SOM_Silhouette-Score.png).
In case the number of cluster chosen is different from the one used to train the SOM, one could read-in the trained SOM, and re-do the analysis with a different number of cluster. To read a trained SOM (and skip the training procedure) simply add the option --SOM and specify the location of the file SOM.Rdata within the output folder:

    PathDetect-SOM.r --folder=COORDS --rep=15 --lig=41-45 --type=dRMSD --cutoff=1.2 --nclus=13 --path_clus=dependent --SOM=SOM/SOM_SOM.Rdata --out=SOM_clus13

The number of cluster (13 in this case) is specified with the option --nclus. All the images will be re-generated in another directory (SOM_clus13 as specified with the option --out).

4. ## **Extract representative conformations**

At the end of the training procedure, each frame of the trajectories is assigned to a neuron of the map. Frames belonging to each neuron can be easily extracted using files within the "Neurons" sub-directory. To extract frames belonging to each neuron simply use the bash script "Extract_Frames.sh". You'll need a trajectory file (xtc or trr) of the simulation, with the same number of frames used to train the SOM. To generate this file simply do within the SIMS folder:[^2] 

    gmx trjcat -f *.xtc -o Full_SIMs.xtc -cat

The simulations will be concatenated in alphabetical order, so make sure that this is the same order of the xvg file in the COORDS folder.
Open the "Extract_Frames.sh" file with a text editor and specify in the first two lines the location (please use absolute path and not relatives) of the concatenated trajectory and a reference file (.gro or .pdb). For example:

    SIM=/home/stefano/work/Tutorial/SIMS/Full_SIMs.xtc
    GRO=/home/stefano/work/Tutorial/SIMS/Protein-lig.gro

Then simply run the bash script to extract frames belongin to each neuron in each corresponding sub-directory:

    bash Extract_Frames.sh

One may also be interested in having a single representative conformation for each neuron. To extract conformation representative of each neuron use the bash script "Extract-Representatives.sh". Open the file with a plain text editor and specify the location of the concatenated trajectory and a reference file as done with the previous file. Run the bash script to extract representative conformations:

    bash Extract-Representatives.sh

This will create a number of .xtc files (e.g. NEURON_0001.xtc) containing a single frame corresponding to the representative conformation for that neuron. You can convert the .xtc files in pdb with the following commands (change the path to the reference gro file in the second line):

    mkdir Representatives
    cp /home/stefano/work/Tutorial/SIMS/Protein-lig.gro Representatives
    mkdir Representatives/XTC
    mv NEURON_*.xtc Representatives/XTC
    cd Representatives
    for FILE in XTC/NEURON_*.xtc; do
    gmx trjconv -f ${FILE} -s Protein-lig.gro -o ${FILE%.*}.pdb <<EOC
    0
    EOC
    done
    mv XTC/*.pdb .



[^1]: Callea, L.; Bonati, L.; Motta, S. Metadynamics-Based Approaches for Modeling the Hypoxia-Inducible Factor 2α Ligand Binding Process. Journal of Chemical Theory and Computation 2021, 18. https://doi.org/10.1021/acs.jctc.1c00114.
[^2]: If you used the option -skip, you'll need to apply the same stride to the trajectory!
