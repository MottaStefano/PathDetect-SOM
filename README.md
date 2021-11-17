# PathDetect-SOM: 
# A tool to analyze pathways from MD simulation on a trained SOM.
A tool based on Self-Organizing Maps (SOM) for the analysis of pathways sampled in molecular dynamics simulations.

1. ## **Requirements**
In order to use the R tool, you will need:

- R (≥ 3.4)
- The following R libraries:
  - kohonen (≥ 3.01)
  - parallel (≥ 3.4)
  - cluster (≥ 2.0.6)
  - igraph (≥ 1.2.4)
- Gromacs

R is available at: [www.r-project.org](http://www.r-project.org)

To install R libraries simply open the R console by typing R on the terminal and run:

    install.packages(kohonen)
    install.packages(parallel)  
    install.packages(cluster)
    install.packages(igraph)

Gromacs is available at: [https://www.gromacs.org/](https://www.gromacs.org)

2. ## **Install PathDetect-SOM**
To install PathDetect-SOM simply copy the script in a directory in your $PATH (e.g. /usr/local/bin in linux) or run from local with ./PathDetect-SOM.r.
Note that in the first line of the code it is specified to run the script with /usr/bin/Rscript, wich is the typical installation folder for Rscript. If Rscript is installed in another folder on your computer, please modify this line.


3. ## **Generate your first Gromacs Coordinate File**
Once you have run MD simulations with your favorite MD engine, you need to extract a set of coordinates that will be used to evaluate the similarity among frames for the SOM training. First of all you need your simulations to be in the gromacs format (.xtc or .trr) and a reference structure file containing the same atoms as in your trajectory files in a format readable by gromacs (.pdb, .gro or .tpr). If you did not use gromacs to generate your trajectory you can use vmd or a tool for trajectory conversion such as mdconvert. Then you must ensure that the periodic boundary conditions are correctly solved in your simulation. In gromacs this will be accomplished with the gmx trjconv command using the option -pbc mol (the best option if you have a single solute molecule in your system but require a tpr file) or -pbc nojump. Always check you output simulation with vmd to be sure that the pbc are correctly solved.

The atoms you choose for PathDetect-SOM must well describe the conformational changes for which you want to track the pathways. Let us assume you followed the unfolding of a protein through multiple replicas[^1]of steered MD simulation (SMD). In this case the conformational change involves the entire protein, so we should choose a set of atoms describing the whole protein. One obvious choice would be to use the set of Cα atoms of the protein. To extract the atom coordinates for frames of your simulations you can use the gromacs gmx traj command:

    gmx traj -f SIM-01.xtc -s Protein.gro -ox coord-01.xvg 

where with the option -f you select the trajectory file of one of your SMD replicas, with option -s you specify the structure file for your system and with option -ox the output coordinate file. Once you have run the command you will be asked for a choice of the atoms you want to extract the coordinates. The choice menu should be something like this:

    Group     0 (         System) has  1782 elements
    Group     1 (        Protein) has  1782 elements
    Group     2 (      Protein-H) has   904 elements
    Group     3 (        C-alpha) has   108 elements
    Group     4 (       Backbone) has   324 elements
    Group     5 (      MainChain) has   433 elements
    Group     6 (   MainChain+Cb) has   537 elements
    Group     7 (    MainChain+H) has   539 elements
    Group     8 (      SideChain) has  1243 elements
    Group     9 (    SideChain-H) has   471 elements

Select a group:

Type “3” and press enter to choose the C-alpha atoms group. If you want to select a group of atoms that is not within the default groups, you can generate an index file with a personalized group of atoms. To do so use the command:

    gmx make_ndx -f Protein.gro -o index.ndx 

Now you can select a group of residues or atoms combining the commands of gmx make\_ndx. Atoms can be selected by name (e.g. a CB), or by number (e.g. a 1-20). You can select a residue (e.g. r 85), you can intersect two existing group selecting only atoms belonging to both groups (e.g. 10 & 3) or combine two groups (e.g. 10 | 11). Pay attention at this selection process because it is easy for beginners to make mistakes.

If your replicas are named with progressive numbers, you can set up a *for* cycle in bash. For example, if you have 20 replicas named SIM-01.xtc, SIM-02.xtc and so on you can do:

    for REP in {01..20}; do
    gmx traj -f SIM-${REP}.xtc -s Protein.gro -ox coord-${REP}.xvg -n index.ndx <<EOC
    3
    EOC
    done

Selecting for all the replicas the group 3 present in the index.ndx file as output group.

4. ## **Prepare your files and run PathDetect-SOM**
Once you have generated the xvg files containing the atom coordinates for the group of your interest, put all this files in a folder.

    mkdir coords/
    mv coord-\*.xvg coords/

PathDetect-SOM will read all the xvg files present in the directory you will indicate, in alphabetical order.

At this point you can run the PathDetect-SOM tool. Running PathDetect-SOM without arguments or with the option --help you obtain the list of arguments accepted by the tool and their meaning. To run the script you must provide at least the folder containing the xvg files and the number of replicas you are analyzing:

    PathDetect-SOM.r --folder=coords --rep=20

In this way the tool assume that all the replicas are of the same length. If you have replicas with different length you can provide the length of the replicas instead of the number of replicas in the option –-rep. In this case specify which would be the first frame of each replica if they would be concatenated. If, for example, you have 10 replicas (1-10) of 1000 frames each and 10 replicas (11-20) of 1500 frames you could provide the option: 

    --rep=1,1001,2001,3001,4001,5001,6001,7001,8001,9001,10001,11501,13001,14501,16001,17501,19001,20501,22001,13501

Other options accepted by PathDetect-SOM tool are:

- --skip to specify a stride in reading the coord files (read 1 line every *n* lines)[^2].
- --dim to specify the dimension of the squared SOM
- --topo to specify the topological shape of the neurons (change the number of neighbors of each neurons.
- --periodic to specify if the SOM would be periodic (neurons on left boundaries neighbors to the ones on the right boundaries and the same for the top and bottom boundaries) or not.
- --dist to specify the distance function used during the SOM calculation
- --tr\_step to specify the number of training cycle for the SOM training
- --mode to specify the type of learning algorithm. The difference between the online and batch SOM algorithms is that the update of the winning unit(s) in the online algorithm is done after each individual object, whereas the batch algorithm does not change the codebook vectors until the whole data set has been presented. pbatch is the parallel version of the batch algorithm.
- --ncores to specify the number of cores to be used during training (if pbatch mode is used)
- --nclust to specify the number of clusters draw on the SOM. Use a guess number the first time, select the optimal number and rerun the command providing as input the trained SOM.
- --clus\_met to specify the clustering method for neuron clustering.
- --dist\_clus to specify the distance to be used for the clustering of neurons.
- --path\_clus to specify the type of pathway clustering. It could be "independent" or "dependend" from time. In the case of time dependent clustering, distances are computed between frames at the same time in the simulations, while with independent clustering distances are computed between frames of a simulation, and the clostest frame of the other simulation.
- --colors to specify a file containing containing a set of 15 hex colors (one per line) to be used for the cluster colors in figures. By default a set of colors is used.
- --out to specify the prefix of all files and output folder.
- --type to specify the type of distance computed. Using RMSD coordinates would be presented to the SOM as the are (simulations must be aligned before the gmx traj command). Using dRMSD the distance matrix between all pair of atoms in the set is computed and this set of distances are used to train the SOM.
- --cont when used together with –-type=dRMSD, allows to consider only distances between atoms making contact in the first frame of the first simulation. This is useful in unfolding process to consider only distances between native contacts.
- --lig In case you are studying a ligand binding/unbinding process with the option –-type=dRMSD, you can use this option to consider only the intermolecular distances. Specify here the ligan atom numbers within the selection (e.g. the selection is composed by 72 atoms of the protein binding site followed by 15 ligand atoms, specify here 73-87).[^3] 
- --cutoff In case you are using the option –-type=dRMSD, distances greater than this value are set at the cutoff value (capping).
- --data A .Rdata object file contained the distance or the coordinates file would be read as input avoiding PathDetect-SOM to re-read the xvg files.
- --SOM A .Rdata object file containing a trained SOM would be read as input. It is useful in case you have already computed the SOM once and you want to simply modify a color, or the number of clusters.
- --seed A value for the random seed. This is usefull for reproducibility. By default a random number will be used.
                                                                            

[^1]: The use of multiple SMD replicas is advisable due to the non-deterministic nature of SMD
[^2]: If you have replicas with different length, adjust also the –rep option accounting for the stride you are applying.
[^3]: This range refers to the atoms within the index, and not to the PDB atom number!                     
                                                                                  
