#!/usr/bin/Rscript

###Parser

## Collect arguments
args <- commandArgs(TRUE)
 
## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat('
      PathDetect-SOM.r
      
      With this command you can train a SOM over data present in a specific folder. This data must be atom coordinates
      in the xvg format extracted with the GROMACS gmx traj command. You can use this script to compare patways sampled
      in multiple replicas of simulations (like Steered MD). In this case you should provide the number of replicas with
      --rep (if the simulations have the same number of frames) or the first frame of each replica considering all
      simulations contatenated.
      
      Mandatory Arguments:
      --folder    - The Input folder for the xvg files. All xvg files present in the folder will be read alphabetically.
      --rep       - The number of replicas performed. If they are not of the same length a comma separated (no space) 
                    string with the frame at wich sepate paths can be supplied. (Considering all simulations 
                    concatenated with the skip) e.g. --rep 1,201,411,621,751.
      
      Optional Arguments:
      
      --skip      = 1             - Read lines with stride when reading coord files
      --dim       = 10            - Dimension of the square SOM in neurons (default 8x8)
      --topo      = hexagonal     - Choose between a hexagonal or rectangular shape of the neuron
      --periodic  = FALSE         - Choose wether the SOM grid is periodic across the boundaries or not (TRUE|FALSE)
      --dist      = "euclidean"   - Distance function to be used for the SOM calculation. Admissable values are: 
                                    "sumofsquares", "euclidean", "manhattan", and "tanimoto".
      --tr_step   = 5000          - Length of SOM training
      --mode      = "pbatch"      - Type of learning algorithm: "online", "batch" or "pbatch".
      --ncores    = -1            - Number of cores to be used in caso of pbatch algorithm. 
                                    By default -1 is all available cores.
      --nclust    = 10            - The number of clusters draw on the SOM
      --clus_met  = complete      - Type of hclust clustering method to be applied on neuron (complete, single, 
                                    average, mcquitty, median, centroid)
      --dist_clus = euclidean     - The distance passsed to hclust for the clustering of neuron. This must be one of:
                                    euclidean, maximum, manhattan, canberra, binary or minkowski.
      --path_clus = independent   - The type of pathway clustering. It could be "independent" or "dependend" from time.
                                    In the case of time dependent clustering, distances are computed between frames at the same
                                    time in the simulations, while with independent clustering distances are computed between frames 
                                    of a simulation, and the clostest frame of the other simulation.
      --colors    = default       - A file containing a set of 15 hex colors (one per line) to be used for the 
                                    cluster colors in figures. By default a set of 13 colors is used.
      --out       = SOM           - Output folder and prefix for the files (if do not exist is created)
      --type      = dRMSD          - Type of distance:
                                    RMSD (sims needs to be pre-aligned) or dRMSD (molecules must be whole)
      --cont      = 0             - If a number (in nm) is given dRMSD, is computed only on distances between atoms forming
                                    contacts in the first frame.
      --lig       = FALSE         - If a range of numbers corresponding to the atoms of ligands within the selection used by 
                                    gmx traj is given (e.g. 181-199), the dRMSD is computed using only intermolecular distances.                     
      --cutoff    = 0             - If a number (in nm) is given, distances greater than this value are set at the cutoff value (capping)
      --data      = no            - If a COORD.Rdata or DIST.Rdata file is supplied here, load it and recompute SOM
      --SOM       = no            - If a SOM.Rdata file is supplied here, load it and only do the analysis and images
      
      Example:
      PathDetect-SOM.r --folder=coords --rep 1,201,411,621,751 --dim=10 --out=SOM_001 --type=dRMSD --colors colors.dat \n\n')
 
  q(save="no")
}



### Functions

#Function to read command line arguments
readArgs <- function(args){
    ## Parse arguments (we expect the form --arg=value)
    parseArgs <- function(x) strsplit(sub("--", "", x), "=")
    argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
    argsL <- as.list(as.character(argsDF$V2))
    names(argsL) <- argsDF$V1

    ## if some mandatory argument is missing:
    if(is.null(argsL$folder) | is.null(argsL$rep) ) {
    cat("Missing Mandatory Argument. Try PathDetect-SOM.r --help\n")
    q(save="no")
    }

    ##Define default values
    if(is.null(argsL$skip)) {
        argsL$skip="1"
    }
    if(is.null(argsL$dim)) {
        argsL$dim="10"
    }
    if(is.null(argsL$topo)) {
        argsL$topo="hexagonal"
    }
    if(is.null(argsL$periodic)) {
        argsL$periodic="FALSE"
    }
    if(is.null(argsL$dist)) {
        argsL$dist="euclidean"
    }
    if(is.null(argsL$tr_step)) {
        argsL$tr_step="5000"
    }
    if(is.null(argsL$mode)) {
        argsL$mode="pbatch"
    }
    if(is.null(argsL$ncores)) {
        argsL$ncores="-1"
    }
    if(is.null(argsL$nclus)) {
        argsL$nclus="10"
    }
    if(is.null(argsL$clus_met)) {
        argsL$clus_met="complete"
    }
    if(is.null(argsL$dist_clus)) {
        argsL$dist_clust="euclidean"
    }
    if(is.null(argsL$path_clus)) {
        argsL$path_clus="independent"
    }
    if(is.null(argsL$colors)) {
        argsL$colors="default"
    }
    if(is.null(argsL$out)) {
        argsL$out="SOM"
    }
    if(is.null(argsL$type)) {
        argsL$type="dRMSD"
    }
    if(is.null(argsL$cont)) {
        argsL$cont="0"
    }
    if(is.null(argsL$lig)) {
        argsL$lig="FALSE"
    }
    if(is.null(argsL$cutoff)) {
        argsL$cutoff="0"
    }
    if(is.null(argsL$graph)) {
        argsL$graph="FALSE"
    }
    if(is.null(argsL$SOM)) {
        argsL$SOM="no"
    }
    if(is.null(argsL$data)) {
        argsL$data="no"
    }
    return(argsL)
}

#Function to compute distance matrix
Calc_Dist_Mat <- function(VEC){
    N_atm <- (length(VEC)-1)/3
    MAT <- as.matrix((dist(t(matrix(VEC[2:length(VEC)], ncol=N_atm)), method='euclidean', upper=TRUE, diag=TRUE)))
    return(MAT)
}

#Function to select only distances between residues making contacts in first frame of the first file (SELE2)
Features_Sele <- function(FILES, SKIP, LIG, CONT){
    L <- system(paste('grep -n "@" ', FILES[1], ' | tail -n 1', sep=''), intern=TRUE)
    L <- as.integer(strsplit(L, split=":")[[1]][1])
    #Read Coord File
    COORDS <- read.table(FILES[1], skip=L)
    SUB_COORDS <- COORDS[seq(1, dim(COORDS)[1], by=SKIP),]
    #How many atoms and frames (skipped)
    N_atm <- (dim(SUB_COORDS)[2]-1)/3
    VEC <- SUB_COORDS[1,]
    D <- Calc_Dist_Mat(VEC)
    DIST_MAT <- array(D, dim=c(N_atm,N_atm,1))
    #Select only lower triangle of matrix
    LOWTRI <- which(lower.tri(D)==TRUE)
    if(LIG != FALSE){
        LIG_a <- strsplit(LIG, "-")[[1]]
        LIG_atoms <- seq(LIG_a[1], LIG_a[2], by=1)
        if(length(which(LIG_atoms %in% c(1:N_atm) ==FALSE)) > 0){
            cat(paste("Ligands atoms: \nLIG_atoms \nare not in the range 1-", N_atm, "\n", sep=''))
            q(save="no")
        }
        PROT_atoms <- which(c(1:N_atm) %in% LIG_atoms ==FALSE)
        DIST_MAT <- DIST_MAT[PROT_atoms,LIG_atoms,1]
        SELE2 <- which(DIST_MAT<CONT)
    } else {
        SELE <- which(DIST_MAT<CONT)
        SELE2 <- SELE[SELE %in% LOWTRI]
    }
    return(SELE2)
}

#Function to compute distances (in case of type=dRMSD)
Compute_Features <- function(FILENAME, SKIP, LIG, SELE2){
    #Read which line is ending intestation
    L <- system(paste('grep -n "@" ', FILENAME, ' | tail -n 1', sep=''), intern=TRUE)
    L <- as.integer(strsplit(L, split=":")[[1]][1])
    #Read Coord File
    COORDS <- read.table(FILENAME, skip=L)
    SUB_COORDS <- COORDS[seq(1, dim(COORDS)[1], by=SKIP),]
    N_atm <- (dim(SUB_COORDS)[2]-1)/3
    #Compute distance matrix for all the frames
    D <- apply(SUB_COORDS, 1, Calc_Dist_Mat)
    if(LIG != FALSE){
        LIG_a <- strsplit(LIG, "-")[[1]]
        LIG_atoms <- seq(LIG_a[1], LIG_a[2], by=1)
        PROT_atoms <- which(c(1:N_atm) %in% LIG_atoms ==FALSE)
        DIST_MAT <- array(D, dim=c(N_atm,N_atm,ncol(D)))
        D <- DIST_MAT[PROT_atoms,LIG_atoms,]
        D <- array(D, dim=c(dim(D)[1]*dim(D)[2], dim(D)[3]))
    }
    cat(sprintf('\r Read file: %s                      ', FILENAME))
    return(D[SELE2,])
}

#Function to read the coordinates (in case of type=RMSD)
Read_coordinates <- function(FILENAME, SKIP){
    #Read which line is ending intestation
    L <- system(paste('grep -n "@" ', FILENAME, ' | tail -n 1', sep=''), intern=TRUE)
    L <- as.integer(strsplit(L, split=":")[[1]][1])
    #Read Coord File
    CCC <- read.table(FILENAME, skip=L)
    COORDS <- rbind(COORDS, CCC[seq(1, dim(CCC)[1], by=SKIP),])
}

#Function to define colorscale
SetColors <- function(COLORS, NCLUS){
    if(COLORS != "default"){
        COL.SCALE <- readLines(COLORS)
        if(length(COL.SCALE < NCLUS)){
            cat(paste("WARNING: Not enought colors in ", COLORS, ", using default\n", sep=''))
            COL.SCALE <- rainbow(NCLUS)
        }
    } else {
        if(NCLUS < 13){
            COL.SCALE <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#a0451f", "#96c3dc", "#a4db77", "#ffff88", "#bea0cc", "#747474", "#f88587", "#fbb25c")
        } else {
            COL.SCALE <- rainbow(NCLUS)
        }
    }
    return(COL.SCALE)
}

#Function to plot circles over a SOM, proportional to a certain property
SOM_circle <- function(FILENAME, SOM, DIM, COL.SCALE, SOM.hc, PROPERTY, TITLE){
    COL <- COL.SCALE[SOM.hc]
    MAGNIF <- 40/DIM
    PROPERTY <- PROPERTY/max(PROPERTY)
    png(FILENAME, width=2200, height=2200)
        par(cex=4)
        plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main=TITLE)
        add.cluster.boundaries(SOM, SOM.hc, lwd=15)
        points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=PROPERTY*MAGNIF*1.2, col="black", xpd=TRUE)
        points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=PROPERTY*MAGNIF, col="white", xpd=TRUE)
    invisible(dev.off())
}

#Function to compute and plot the silhouette profiles
Do_Silhouette <- function(PREFIX, SOM, DIST_CLUS, CLUS_METHOD){
    SIL <- NULL
    for(i in 2:50){
        png(paste(PREFIX, "_Silhouette-", i, ".png", sep=''), width=2200, height=2200)
            par(cex=4, lwd=3)
            SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), i)
            sil = silhouette(SOM.hc, dist(SOM$codes[[1]]))
            plot(sil)
            abline(v=mean(sil[,3]), col='red', lty=3, lwd=2)
        invisible(dev.off())
        SIL <- c(SIL, mean(sil[,3]))
    }
    return(SIL)
}

#Function to compute the weighted mean (by population) of the vectors belonging to each clusters
clust.centroid = function(i, dat, clusters) {
    ind = (clusters == i)
    if(sum(ind)>1){
        POP <- NULL
        for(NEURON in 1:nrow(SOM$grid$pts)){
            POP <- c(POP, length(which(SOM$unit.classif==NEURON)))
        }
         apply(dat[ind,], 2, weighted.mean, w=POP[ind])
    } else {
        dat[ind,]
    }
}

#Function to compute the distance between two vectors
COMPUTE_DISTANCE <- function(V1, V2){
    return(as.numeric(dist(rbind(V1, V2))))
}

#Function to create a vector of two letters identifiers
TwoLetters <- function(LETTERS=LETTERS){
    LET <- NULL
    for(i in LETTERS){
        for(j in LETTERS){
        LET <- c(LET, paste(i, j, sep=''))
        }
    }
    return(LET)
}

#Function to create the cluster names
Define_ClusName <- function(NCLUS){
    if(NCLUS < length(LETTERS)){
        CLUS_NAME <- paste("Cluster", LETTERS, sep="_")[1:NCLUS]
    } else{
#         LET <- TwoLetters(LETTERS)
        CLUS_NAME <- paste("Cluster", LET, sep="_")[1:NCLUS]
    }
    return(CLUS_NAME)
}

#Function to select the neuron representative of cluster CL
Select_representative <- function(CENTROID, SOM, SOM.hc, CL){
    FR <- which(SOM.hc==CL)
    DistCentroid <- apply(SOM$codes[[1]], 1, COMPUTE_DISTANCE, V2=CENTROID[,CL])
    ReprNeuron <- which(DistCentroid==min(DistCentroid[FR]))
    return(ReprNeuron)
}

#Function to compute the representative frame for each neuron
NeuronRepres <- function(SOM){
    REPRESENTATIVE <- NULL
    for(i in 1:nrow(SOM$grid$pts)){
        SET <- which(SOM$unit.classif==i)
        if(length(SET)==0){
            REPRESENTATIVE <- c(REPRESENTATIVE, 0)
        } else{
            REPRESENTATIVE <- c(REPRESENTATIVE, SET[which(SOM$distances[SET]==min(SOM$distances[SET]))])
        }
    }
    return(REPRESENTATIVE)
}

#Function to write a series of frames (NUMs) to a file in the ndx format
WriteNDX <- function(FILENAME, TITLE, NUMs){
    write(TITLE, file=FILENAME)
    for (j in 0:(as.integer(length(NUMs)/15)-1)){
        L <- NUMs[(j*15)+1]
        for(k in ((j*15)+2):((j+1)*15)){
            L <- paste(L, NUMs[k], sep=' ')
        }
        write(L, file=FILENAME, append=TRUE)
    }
    L <- NULL
    if((j+1)*15<length(NUMs)){
        for(k in NUMs[c((((j+1)*15)+1):length(NUMs))]){
            L <- paste(L, k, sep=' ')
        }
    }
    write(L, file=FILENAME, append=TRUE)
}

#Function to plot a SOM with neuron numbering
plot_SOM_Number <- function(SOM, SOM.hc, COL){
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main="SOM Numbering")
    add.cluster.boundaries(SOM, SOM.hc, lwd=15)
    X <- NULL
    Y <- NULL
    for(i in c(1:nrow(SOM$grid$pts))){
        X <- SOM$grid$pts[i,1]
        Y <- SOM$grid$pts[i,2]
        text(x=X, y=Y, labels=i, cex=(13/SOM$grid$xdim), xpd=TRUE)
    }
}

#Function to apply a legend of clusters to a SOM map image
Legend_Cluster <- function(NCLUS, COL.SCALE){
    if(NCLUS < length(LETTERS)){
        LEG_LAB <- paste("Cluster ", LETTERS, sep=" ")[1:NCLUS]
        CX <- 1.24
    } else{
        LET <- TwoLetters(LETTERS)
        LEG_LAB <- paste("Cluster ", LET, sep=" ")[1:NCLUS]
        CX <- 1.24
    }
    NR=ceiling(NCLUS/4)
    MyOrder = matrix(1:(NR*4), nrow=NR, ncol=4, byrow=T)
    MyBorders = rep("black", NR*4)
    MyBorders[MyOrder > NCLUS] <- NA
    legend(x=0.4, y=0, legend=LEG_LAB[MyOrder], fill=COL.SCALE[1:NCLUS][MyOrder], ncol=4, xpd=TRUE, cex=CX, bty="n", border=MyBorders) #or pch=22
}

#Function to draw pathways over the SOM
Trace_Path <- function(SOM, N, STARTS, ENDS, STRIDE=1){
    X <- NULL
    Y <- NULL
    BWR <- colorRampPalette(c("blue", "white", "red"))
    for(i in seq(STARTS[N], ENDS[N], by=STRIDE)){
        u <- SOM$unit.classif[i]
        X <- c(X, SOM$grid$pts[u,1])
        Y <- c(Y, SOM$grid$pts[u,2])
    }
    points(X,Y, pch=16, cex=30/SOM$grid$xdim)
    points(X,Y, pch=16, col=BWR(length(X)), cex=24/SOM$grid$xdim)
    lines(X,Y, pch=16, lwd=5)
}

#Function to compute the distance between two paths in a time-dependent mode
dist_path_TimeDep <- function(A, B, GRID){
    #A and B are two vectors of the same length containing the path through neurons while grid is the SOM$grid$pts
    D <- NULL
    for(i in 1:length(A)){
        D <- c(D, dist(rbind(GRID[A[i],], GRID[B[i],]), method="euclidean"))
    }
    return(sum(D)/length(A))
}

#Function to compute the distance between two paths in a time-independent mode
dist_path_TimeIndep <- function(A, B, GRID){
    #A and B are two vectors of the same length containing the path through neurons while grid is the SOM$grid$pts
    D1 <- NULL
    for(i in 1:length(A)){
        D1 <- c(D1, min(as.matrix(dist(rbind(GRID[A[i],], GRID[unique(B),]), method="euclidean", upper=TRUE, diag=TRUE))[1,-1]))
    }
    D2 <- NULL
    for(i in 1:length(B)){
        D2 <- c(D2, min(as.matrix(dist(rbind(GRID[B[i],], GRID[unique(A),]), method="euclidean", upper=TRUE, diag=TRUE))[1,-1]))
    }
    D1 <- sum(D1)/length(A)
    D2 <- sum(D2)/length(B)
    return(max(D1,D2))
}

#Compute the transition matrix starting from a vector of subsequent classifications
Compute_Transition_Matrix <- function(CLASSIF){
    #Compute the probability of passing from neuron i to neuron j
    TRANS <- matrix(0, ncol=max(CLASSIF), nrow=max(CLASSIF))
    for(i in 1:max(CLASSIF)){
        #Total number of frame assigned to the neuron i
        TOTAL <- length(which(CLASSIF==i))
        #Neurons to which the neuron i has evolved to
        PASSAGE <- CLASSIF[which(CLASSIF==i)+1]
        if(TOTAL > 0){
            for(j in 1:max(CLASSIF)){
                NNN <- length(which(PASSAGE==j))
                    TRANS[i,j] <- NNN
            }
        }
    }
    colnames(TRANS) <- paste("N_", seq(1:max(CLASSIF)), sep='')
    rownames(TRANS) <- paste("N_", seq(1:max(CLASSIF)), sep='')
    return(TRANS)
}

#Function to convert a transition matrix to a network
Matrix2Network <- function(TRANS){
    #Create transition network
    d <- NULL
    for(i in 1:dim(TRANS)[1]){
        for(j in 1:dim(TRANS)[2]){
            if(TRANS[i,j]>0){
                d <- rbind(d, c(i, j, TRANS[i,j]))
            }
        }
    }
    return(d)
}

#Function to remove the diagonal from a network
Network_noDiagonal <- function(d){
    #Remove elements on the diagonal
    D <- NULL
    for(i in 1:nrow(d)){
        if(d[i,1]!=d[i,2]){
            D <- rbind(D, d[i,])
        }
    }
    return(D)
}

#Function to convert a network dataframe to an igraph graph
Network2Graph <- function(D, SOM, SOM.hc, COL.SCALE){
    POPULATION <- NULL
    N_NEUR <- nrow(SOM$codes[[1]])
    for(i in 1:N_NEUR){
        POPULATION <- c(POPULATION, length(which(SOM$unit.classif==i)))
    }
    #Create dataframe for nodes
    nodes <- cbind(c(1:nrow(SOM$grid$pts)), SOM.hc, POPULATION)
    colnames(nodes) <- c("node", "CLUSTER", "POPULATION")
    #Create Igraph Network
    colnames(D) <- c("V1","V2","weight")
    net <- graph_from_data_frame(d=D, vertices=nodes, directed=T)
    #Set some properties of the graph
    V(net)$color <- COL.SCALE[V(net)$CLUSTER]
    V(net)$size <- log(V(net)$POPULATION)*1.2
    Isolated = which(degree(net)==0)
    V(net)$size[Isolated] <- 0
    E(net)$width <- (log(D[,3])/max(log(D[,3])))*10
    return(net)
}



#################################################################################################################
#################################################################################################################
#########################################              Main             #########################################
#################################################################################################################
#################################################################################################################

############################### Read Command line and assign values to variables ################################

argsL <- readArgs(args)
# Values from args to variables
FOLDER      <- argsL$folder
REP         <- strsplit(argsL$rep, ",")[[1]]
SKIP        <- as.integer(argsL$skip)
DIM         <- as.integer(argsL$dim)
TOPO        <- argsL$topo
PERIODIC    <- argsL$periodic
DISTANCE    <- argsL$dist
RLEN        <- as.integer(argsL$tr_step)
MODE        <- argsL$mode
NCORES      <- as.integer(argsL$ncores)
NCLUS       <- as.integer(argsL$nclus)
CLUS_METHOD <- argsL$clus_met
DIST_CLUS   <- argsL$dist_clus
PATH_CLUS   <- argsL$path_clus
COLORS      <- argsL$colors
OUT         <- argsL$out
TYPE        <- argsL$type
CONT        <- as.numeric(argsL$cont)
LIG         <- argsL$lig
CUT         <- as.numeric(argsL$cut)
GRAPH       <- argsL$graph
SOM_file    <- argsL$SOM
DATA        <- argsL$data
if(length(REP)==1){
    NREP <- as.integer(REP)
} else {
    NREP <- length(REP)
}
if(CONT==0){
    CONT <- "Inf"
}
#Write Params to file
dir.create(OUT, showWarnings=FALSE)
PARAMS <- cbind(c("--folder     = ", "--rep        = ", "--skip       = ",
                  "--dim        = ", "--topo       = ", "--periodic   = ",
                  "--dist       = ", "--tr_step    = ", "--mode       = ",
                  "--ncores     = ", "--nclust     = ", "--clus_met   = ",
                  "--dist_clus  = ", "--path_clus  = ", "--colors     = ",
                  "--out        = ", "--type       = ", "--cont       = ",
                  "--lig        = ", "--cutoff     = ", "--graph      = ",
                  "--SOM        = ", "--data       = "),
                c(FOLDER , paste(REP, collapse=',') , SKIP , DIM , TOPO ,
                  PERIODIC , DISTANCE , RLEN , MODE , NCORES , NCLUS ,
                  CLUS_METHOD, DIST_CLUS , PATH_CLUS , COLORS , OUT ,
                  TYPE , CONT , LIG , CUT , GRAPH , SOM_file , DATA))
write.table(PARAMS, file=paste(OUT, "/", OUT, "_PARAMETERS.dat", sep=''), row.names=FALSE, col.names=FALSE, quote=FALSE)


################################# FEATURES CALCULATION ##################################
if(SOM_file=="no"){
    if(DATA=="no"){
        FILES <- system(paste("ls ", FOLDER, "/*.xvg", sep=''), intern=TRUE)
        if(TYPE=="dRMSD"){
            #Define the distances that should be included as features
            cat(sprintf(' Selecting Features...\n'))
            FSele <- Features_Sele(FILES, SKIP, LIG, CONT)
            DIST <- NULL
            #Compute distances for all the files
            for(N in 1:length(FILES)){
                cat(sprintf('\r Parsing file: %s', FILES[N]))
                DIST <- cbind(DIST, Compute_Features(FILES[N], SKIP, LIG, FSele))
            }
            cat(sprintf('\n Done'))
            #Apply capping
            if(CUT>0){
                DIST[which(DIST > CUT)] <- CUT
            }
            #Rotate the matrix
            DIST <- t(DIST)
            cat(sprintf('\n Saving DIST Matrix..'))
            save(DIST, file=paste(OUT, "DIST.Robj", sep=''))
        }
        if(TYPE=="RMSD"){
            #Read and store atom coordinates
            COORDS <- NULL
            for(N in 1:length(FILES)){
                cat(sprintf('\r Parsing file: %s', FILES[N]))
                Read_coordinates(FILES[N], SKIP)
            }
            cat(sprintf('\n Saving COORD Matrix..'))
            save(COORDS, file=paste(OUT, "COORDS.Robj", sep=''))
        }
        cat(sprintf('\n Done'))
    } else{
        load(DATA)
    }
}

######################################### SOM TRAINING #############################################
library(parallel)
library(kohonen)

if(SOM_file=="no"){
    cat(sprintf('\n Training SOM...'))
    if(TYPE=="RMSD"){
        INP_DATA <- as.matrix(COORDS[,2:ncol(COORDS)])
    }
    if(TYPE=="dRMSD"){
        INP_DATA <- as.matrix(DIST)
    }
    if(PERIODIC=="FALSE"){
        SOM <- som(as.matrix(INP_DATA), grid = somgrid(DIM, DIM, TOPO, neighbourhood.fct="gaussian"), dist.fcts=DISTANCE, rlen=RLEN, mode='pbatch')
    } else {
        SOM <- som(as.matrix(INP_DATA), grid = somgrid(DIM, DIM, TOPO, neighbourhood.fct="gaussian", toroidal=TRUE), dist.fcts=DISTANCE, rlen=RLEN, mode='pbatch')
    }
    cat(sprintf('\n Done'))
    cat(sprintf('\n Saving SOM...'))
    save(SOM, file=paste(OUT, "/", OUT, "_SOM.Rdata", sep=''))
    cat(sprintf('\n Done'))
} else{
    cat(sprintf('\n Reading SOM...'))
    load(SOM_file)
}

#Clustering of neurons
SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)


############################################ BASIC ANALYSIS ##############################################

library(cluster)

### TRAINING
COL.SCALE <- SetColors(COLORS, NCLUS)
cat("\n Creating images for training...\n")
FOLDER <- paste(OUT, "/", "Training", sep='')
dir.create(FOLDER, showWarnings=FALSE)
png(paste(FOLDER, "/", OUT, "_Convergency.png", sep=''), width=2200, height=2200)
par(cex=3.5, lwd=4.5)
plot(SOM, type = "changes", main="Convergency", cex=2.5, lwd=4)
invisible(dev.off())
png(paste(FOLDER, "/", OUT, "_heatmap.png", sep=''), width=2200, height=2200)
par(cex=4)
SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)
plot(SOM, type = 'dist.neighbours', heatkey = TRUE, shape='straight', main="Distance Heatmap")
add.cluster.boundaries(SOM, SOM.hc, lwd=8)
invisible(dev.off())
POP <- NULL
for(NEURON in 1:nrow(SOM$grid$pts)){
    POP <- c(POP, length(which(SOM$unit.classif==NEURON)))
}
#Plot circles proportional to the population of neurons
SOM_circle(FILENAME=paste(FOLDER, "/", OUT, "_population.png", sep=''),
            SOM=SOM, DIM=DIM, COL.SCALE=COL.SCALE, SOM.hc=SOM.hc, PROPERTY=POP, TITLE="Population")
#Compute the population on a logarithm scale
logPOP <- log(POP/length(SOM$unit.classif))
M <- min(logPOP[is.finite(logPOP)])
logPOP <- logPOP-(M-1)
logPOP[which(logPOP=="-Inf")] <- 0
SOM_circle(FILENAME=paste(FOLDER, "/", OUT, "_log-population.png", sep=''),
            SOM=SOM, DIM=DIM, COL.SCALE=COL.SCALE, SOM.hc=SOM.hc, PROPERTY=logPOP, TITLE="log-Population")
cat(" Done\n")



### POPULATION PER-REPLICA
cat(" Creating Per-replica Population Images...\n")
FOLDER <- paste(OUT, "/Population", sep='')
dir.create(FOLDER, showWarnings=FALSE)

Tot_Frame <- round(length(SOM$unit.classif))
if(length(REP)==1){
    STARTS <- seq(1, Tot_Frame, by=Tot_Frame/NREP)
    ENDS <- seq((Tot_Frame/NREP), Tot_Frame, by=Tot_Frame/NREP)
} else {
    REP <- as.numeric(REP)
    STARTS <- REP
    ENDS <- c(STARTS[2:length(STARTS)]-1, Tot_Frame)
}

for(SIM in sprintf("%04d", c(1:NREP))){
    N <- as.numeric(SIM)
    COL <- COL.SCALE[SOM.hc]

    POP <- NULL
    for(NEURON in 1:nrow(SOM$grid$pts)){
        POP <- c(POP, length(which(SOM$unit.classif[STARTS[N]:ENDS[N]]==NEURON)))
    }
    SOM_circle(paste(FOLDER, "/", OUT, "_Rep_", SIM, "_population.png", sep=''),
                SOM=SOM, DIM=DIM, COL.SCALE=COL.SCALE, SOM.hc=SOM.hc, PROPERTY=POP, TITLE=sprintf("Replica %04d", N))
    logPOP <- log(POP/length(SOM$unit.classif[STARTS[N]:ENDS[N]]))
    M <- min(logPOP[is.finite(logPOP)])
    logPOP <- logPOP-(M-1)
    logPOP[which(logPOP=="-Inf")] <- 0
    dir.create(paste(OUT, "/Population/Log", sep=''), showWarnings=FALSE)
    SOM_circle(paste(FOLDER, "/Log/", OUT, "_Rep_", SIM, "_population.png", sep=''),
                SOM=SOM, DIM=DIM, COL.SCALE=COL.SCALE, SOM.hc=SOM.hc, PROPERTY=logPOP, TITLE=sprintf("Replica %04d", N))
}
cat(" Done\n")



### Clusters
cat(" Creating images at different cluster levels...\n")
FOLDER <- paste(OUT, "/", "Clusters", sep='')
dir.create(FOLDER, showWarnings=FALSE)
#Do Images for clusters in range 1-15
for(i in 2:15){
    SOM.hc_i <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), i)
    png(paste(FOLDER, "/", OUT, "_CLUS-", i, ".png", sep=''), width=2200, height=2200)
    par(cex=4)
    if(length(COL.SCALE) < i){
        COL <- rainbow(i)
    } else{
        COL <- COL.SCALE
    }
    COL <- COL[SOM.hc_i]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight')
    add.cluster.boundaries(SOM, SOM.hc_i, lwd=10)
    invisible(dev.off())
}
cat(" Done\n")
cat(" Creating Silhouette Images...\n")
FOLDER <- paste(OUT, "/Clusters/Silhouette/", sep='')
dir.create(FOLDER, showWarnings=FALSE)
SIL <- Do_Silhouette(paste(FOLDER, "/", OUT, sep=''), SOM, DIST_CLUS, CLUS_METHOD)
#Plot the silhouette score profile
png(paste(FOLDER, "/", OUT, "_Silhouette-Score.png", sep=''), width=2300, height=1700)
    par(cex=4, lwd=2)
    plot(c(2:50), SIL, type='b', pch=19, lwd=2, xlab="Number of clusters", ylab='Average silhouettes', axes = FALSE)
    axis(side = 1, at = seq(0, 50, by=5), lwd=2)
    axis(side = 2, lwd=2)
    box()
invisible(dev.off())
cat(" Done\n")



### Definition of CLUSTER REPRESENTATIVES
cat(" Computing Cluster Representatives...\n")
FOLDER <- paste(OUT, "/", "Clusters", sep='')
CLUS_NAME <- Define_ClusName(NCLUS)
#Compute the center of each cluster (weighted mean of vectors)
CENTROID <- sapply(unique(SOM.hc), clust.centroid, SOM$codes[[1]], SOM.hc)
#Select representative neuron (neuron with vector closer to the average vector)
invisible(suppressWarnings(file.remove(paste(FOLDER, "/Cluster_Representatives.dat", sep=''))))
invisible(file.create(paste(FOLDER, "/Cluster_Representatives.dat", sep='')))
ReprNeurons <- NULL
for(i in 1:ncol(CENTROID)){
    ReprNeurons <- c(ReprNeurons, Select_representative(CENTROID, SOM, SOM.hc, i))
    write(paste(CLUS_NAME[i], ": Neuron ",  ReprNeurons[i], sep=''), file=paste(FOLDER, "/Cluster_Representatives.dat", sep=''), append=TRUE)
}
cat(" Done\n")

### Definition of NEURONS REPRESENTATIVE FRAME
cat(" Generating file for the Extraction of Neuron Centroids...\n")
FOLDER <- paste(OUT, "/Neurons", sep='')
dir.create(FOLDER, showWarnings=FALSE)
#Compute the representative frame for each neuron
REPRESENTATIVE <- NeuronRepres(SOM)
#Write files for the extraction of neurons representative frame
write("#!/bin/bash", file=paste(FOLDER, "/Extract_Frames.sh", sep=''))
write("SIM=", file=paste(FOLDER, "/Extract_Frames.sh", sep=''), append=TRUE)
write("GRO=\n", file=paste(FOLDER, "/Extract_Frames.sh", sep=''), append=TRUE)
invisible(suppressWarnings(file.remove(paste(FOLDER, "/REPRESENTATIVES.ndx", sep=''))))
invisible(file.create(paste(FOLDER, "/REPRESENTATIVES.ndx", sep='')))
for(i in 1:length(REPRESENTATIVE)){
    if(REPRESENTATIVE[i]>0){
        write(sprintf("[ NEURON_%04d ]", i), file=paste(FOLDER, "/REPRESENTATIVES.ndx", sep=''), append=TRUE)
        write(REPRESENTATIVE[i], file=paste(FOLDER, "/REPRESENTATIVES.ndx", sep=''), append=TRUE)
    }
}
write("#!/bin/bash", file=paste(FOLDER, "/Extract-Representatives.sh", sep=''))
write("SIM=", file=paste(FOLDER, "/Extract-Representatives.sh", sep=''), append=TRUE)
write("GRO=\n", file=paste(FOLDER, "/Extract-Representatives.sh", sep=''), append=TRUE)
write(paste("gmx trjconv -f $SIM -s $GRO -sub REPRESENTATIVES.ndx -o REPRESENTATIVES.xtc <<EOC\n0\nEOC", sep=''), file=paste(FOLDER, "/Extract-Representatives.sh", sep=''), append=TRUE)

### Creates Files for the Extraction of FRAMES BELONGING TO EACH NEURON
for(i in 1:nrow(SOM$grid$pts)){
    NUMs <- which(SOM$unit.classif==i)
    if(dir.exists(paste(FOLDER, sprintf("/Neuron_%04d", i), sep=''))==FALSE){
        dir.create(paste(FOLDER, sprintf("/Neuron_%04d", i), sep=''))
    }
    WriteNDX(FILENAME=paste(FOLDER, sprintf("/Neuron_%04d/NEURON.ndx", i), sep=''),
             TITLE=sprintf("[ NEURON_%04d ]", i), NUMs=NUMs)
    write(paste("cd ", sprintf("Neuron_%04d/\n", i), "gmx trjconv -f $SIM -s $GRO -sub NEURON.ndx -o NEURON.xtc <<EOC\n0\nEOC\ncd ..\n",     sep=''), file=paste(FOLDER, "/Extract_Frames.sh", sep=''), append=TRUE)
}

### SOM with Neuron Numbering and cluster legend
H <- 2300+145*as.integer(NCLUS/4)
png(paste(FOLDER, "/", OUT, "_SOM-Numbering.png", sep=''), width=2200, height=H)
    par(cex=4)
    plot_SOM_Number(SOM, SOM.hc, COL.SCALE[SOM.hc])
    Legend_Cluster(NCLUS, COL.SCALE)
invisible(dev.off())

cat(" Done\n")




############################################ PATHWAYS RECONSTRUCTION ##############################################

### Per-replica path images
cat(" Creating Per-replica Paths Images...\n")
FOLDER <- paste(OUT, "/Paths", sep='')
dir.create(FOLDER, showWarnings=FALSE)
#Stride value to depict the path
STRIDE <- 1

for(SIM in sprintf("%04d", c(1:NREP))){
    N <- as.numeric(SIM)
    png(paste(FOLDER, "/", OUT, "_Rep_", SIM, "_Path.png", sep=''), width=2000, height=2000)
    par(cex=4)
    COL <- COL.SCALE[SOM.hc]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main=sprintf("Replica %04d", N))
    add.cluster.boundaries(SOM, SOM.hc, lwd=15)
    Trace_Path(SOM, N, STARTS, ENDS, STRIDE=1)
    invisible(dev.off())
}
cat(" Done\n")

### CLUSTERING of pathways
if(NREP > 2){
    if(length(REP)==1){
    steer_PATHS <- matrix(SOM$unit.classif, ncol=NREP)
    } else {
    REPLICAS <- c(REP, length(SOM$unit.classif))
    L <- NULL
    for(i in 2:length(REPLICAS)){
        L <- c(L, REPLICAS[i]-REPLICAS[i-1])
    }
    steer_PATHS <- matrix(NA, nrow=max(L), ncol=length(REPLICAS))
    for(i in 2:length(REPLICAS)){
        steer_PATHS[c(1:L[i-1]), (i-1)] <- SOM$unit.classif[c(REPLICAS[i-1]:(REPLICAS[i]-1))]
    }
    }
    #Compute the distance matrix
    MAT <- matrix(0, ncol=NREP, nrow=NREP)
    if(PATH_CLUS=="dependent"){
        if(length(REP)==1){
            for(i in 1:NREP){
                for(j in 1:NREP){
                    MAT[i,j] <- dist_path_TimeDep(na.omit(steer_PATHS[,i]), na.omit(steer_PATHS[,j]), SOM$grid$pts)
                }
            }
        } else{
            cat(" WARNING: Specified replicas of different length, using time-independent clustering...\n")
            PATH_CLUS="independent"
        }
    }
    if(PATH_CLUS=="independent"){
            for(i in 1:NREP){
                for(j in 1:NREP){
                    MAT[i,j] <- dist_path_TimeIndep(na.omit(steer_PATHS[,i]), na.omit(steer_PATHS[,j]), SOM$grid$pts)
                }
            }
    }
    #Plot dendrograms
    W=1000+50*NREP
    path.clust <- hclust(as.dist(MAT), method="complete")
    png(paste(FOLDER, "/", OUT, "_PATH_Clustering_dendrogram-complete.png", sep=''), width=W, height=2000)
    par(cex=5, lwd=5)
    plot(path.clust, xlab="", main="Clustering dendrogram complete linkage")
    invisible(dev.off())
    path.clust <- hclust(as.dist(MAT), method="average")
    png(paste(FOLDER, "/", OUT, "_PATH_Clustering_dendrogram-average.png", sep=''), width=W, height=2000)
    par(cex=5, lwd=5)
    plot(path.clust, xlab="", main="Clustering dendrogram average linkage")
    invisible(dev.off())
    path.clust <- hclust(as.dist(MAT), method="single")
    png(paste(FOLDER, "/", OUT, "_PATH_Clustering_dendrogram-single.png", sep=''), width=W, height=2000)
    par(cex=5, lwd=5)
    plot(path.clust, xlab="", main="Clustering dendrogram single linkage")
    invisible(dev.off())
    cat(" Done\n")
}



############################################ NETWORK BUILDING ############################################
suppressMessages(library(igraph))
cat(" Generating Network...\n")
FOLDER <- paste(OUT, "/Transition_Matrix", sep='')
dir.create(FOLDER, showWarnings=FALSE)

#Create vector containing the classification (Neuron) assigned to each Frame
CLASSIF <- SOM$unit.classif
#Replace the first frame of each replica with 0 to avoid the passage from last frame to the first
CLASSIF[STARTS] <- 0

TRANS <- Compute_Transition_Matrix(CLASSIF)
write.table(TRANS, paste(FOLDER, "/", OUT, "_Transition_Matrix.dat", sep=''))
d <- Matrix2Network(TRANS)
write.table(d, paste(FOLDER, "/", OUT, "_Transition_Network.dat", sep=''), col.names=FALSE, row.names=FALSE)
# write_graph(d, file=paste(FOLDER, "/", OUT, "_Transition_Network.xml", sep=''), format="graphml")
D <- Network_noDiagonal(d)
write.table(D, paste(FOLDER, "/", OUT, "_Transition_Network-noDiagonal.dat", sep=''), col.names=FALSE, row.names=FALSE)
net <- Network2Graph(D, SOM, SOM.hc, COL.SCALE)
write_graph(net, file=paste(FOLDER, "/", OUT, "_Transition_Network-noDiagonal.xml", sep=''), format="graphml")
cat(" Done\n")
############################################ SAVE PROJECT ############################################
#Save the project
cat(" Saving Project...\n")
save.image(paste(OUT, "/", OUT, "_SOM-Project.Rdata", sep=''))
cat(" Done\n")
