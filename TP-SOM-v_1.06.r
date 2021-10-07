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
      SOM_Path.r
      
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
      --dim       = 8             - Dimension of the square SOM in neurons (default 8x8)
      --topo      = hexagonal     - Choose between a hexagonal or rectangular shape of the neuron
      --periodic  = FALSE         - Choose wether the SOM grid is periodic across the boundaries or not (TRUE|FALSE)
      --dist      = "euclidean"   - Distance function to be used for the SOM calculation. Admissable values are: 
                                    "sumofsquares", "euclidean", "manhattan", and "tanimoto".
      --tr_step   = 5000          - Length of SOM training
      --mode      = "pbatch"      - Type of learning algorithm: "online", "batch" or "pbatch".
      --ncores    = -1            - Number of cores to be used in caso of pbatch algorithm. 
                                    By default -1 is all available cores.
      --nclust    = 6             - The number of clusters draw on the SOM
      --clus_met  = complete      - Type of hclust clustering method to be applied on neuron (complete, single, 
                                    average, mcquitty, median, centroid)
      --dist_clus = euclidean     - The distance passsed to hclust for the clustering of neuron. This must be one of:
                                    euclidean, maximum, manhattan, canberra, binary or minkowski.
      --path_clus = independent   - The type of pathway clustering. It could be "independent" or "dependend" from time.
                                    In the case of time dependent clustering, distances are computed between frames at the same
                                    time in the simulations, while with independent clustering distances are computed between frames 
                                    of a simulation, and the clostest frame of the second simulation.
      --colors    = default       - A file containing a set of 15 hex colors (one per line) to be used for the 
                                    cluster colors in figures. By default a rainbow function is used.
      --out       = SOM           - Output folder and prefix for the files (if do not exist is created)
      --type      = RMSD          - Type of distance: 
                                    RMSD (sims needs to be pre-aligned) or dRMSD (molecules must be whole)
      --cont      = 0             - If a number (in nm) is given dRMSD, is computed only on distances between atoms forming
                                    contacts in the first frame.
      --lig       = FALSE         - If a range of numbers corresponding to the atoms of ligands within the selection used by 
                                    gmx traj is given (e.g. 181-199), the dRMSD is computed using only intermolecular distances.                     
      --cutoff    = 0             - If a number (in nm) is given, distances greater than this value are set at the cutoff value
      --graph     = FALSE         - Represent paths as graph (require igraph)
      --SOM       = no            - If a SOM.Rdata file is supplied here, load it and only do the images
      --data      = no            - If a COORD.Rdata or DIST. Rdata file is supplied here, load it and recompute SOM
      
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
        argsL$dim="8"
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
        argsL$nclus="6"
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
        argsL$type="RMSD"
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
    cat(sprintf('\r Read file: %s', FILENAME))
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


############################### CALCULATION OF COORDS or INTERNAL DISTANCES ################################
if(SOM_file=="no"){
    if(DATA=="no"){
        FILES <- system(paste("ls ", FOLDER, "/*.xvg", sep=''), intern=TRUE)
        if(TYPE=="dRMSD"){
            #Define the distances that should be included as features
            AtomSele <- Features_Sele(FILES, SKIP, LIG, CONT)
            DIST <- NULL
            #Compute distances for all the files
            for(N in 1:length(FILES)){
                cat(sprintf('\r Parsing file: %s', FILES[N]))
                DIST <- cbind(DIST, Compute_Features(FILES[N], SKIP, LIG, SELE2))
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

######################################### CALCULATION of SOM #############################################
library(parallel)
library(kohonen)

dir.create(OUT, showWarnings=FALSE)
if(SOM_file=="no"){
    cat(sprintf('\n Training SOM..'))
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
    cat(sprintf('\n Saving SOM..'))
    save(SOM, file=paste(OUT, "/", OUT, "_SOM.Rdata", sep=''))
    cat(sprintf('\n Done'))
} else{
    load(SOM_file)
}




############################################ DO IMAGES ##############################################

library(cluster)
cat("Creating images for training\n")

### TRAINING
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

#Magnif is the aplitude of points draw to represent the population of each neuron.
MAGNIF <- 40/DIM
SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)
png(paste(FOLDER, "/", OUT, "_log-population.png", sep=''), width=2200, height=2200)
    par(cex=4)
    COL <- COL.SCALE[SOM.hc]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main="log-Population")
    add.cluster.boundaries(SOM, SOM.hc, lwd=15)
    POP <- NULL
    for(NEURON in 1:nrow(SOM$grid$pts)){
        POP <- c(POP, length(which(SOM$unit.classif==NEURON)))
    }
    POP <- POP/length(SOM$unit.classif)
    POP <- log(POP)
    M <- min(POP[is.finite(POP)])
    POP <- POP-(M-1)
    POP[which(POP=="-Inf")] <- 0
    POP <- POP/max(POP)
    points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=POP*MAGNIF*1.2, col="black", xpd=TRUE)
    points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=POP*MAGNIF, col="white", xpd=TRUE)
invisible(dev.off())

#Magnif is the aplitude of points draw to represent the population of each neuron.
SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)
png(paste(FOLDER, "/", OUT, "-population.png", sep=''), width=2200, height=2200)
    par(cex=4)
    COL <- COL.SCALE[SOM.hc]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main="Population")
    add.cluster.boundaries(SOM, SOM.hc, lwd=15)
    POP <- NULL
    for(NEURON in 1:nrow(SOM$grid$pts)){
        POP <- c(POP, length(which(SOM$unit.classif==NEURON)))
    }
    POP <- POP/max(POP)
    points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=POP*MAGNIF*1.2, col="black", xpd=TRUE)
    points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=POP*MAGNIF, col="white", xpd=TRUE)
invisible(dev.off())

cat("Done\n")



### Clusters
cat("Creating images at different cluster levels...\n")
FOLDER <- paste(OUT, "/", "Cluster_Number", sep='')
dir.create(FOLDER, showWarnings=FALSE)
#Do Images for clusters in range 1-15
for(i in 2:15){
    png(paste(FOLDER, "/", OUT, "_CLUS-", i, ".png", sep=''), width=2200, height=2200)
    par(cex=4)
    if(length(COL.SCALE) < i){
        COL <- rainbow(i)
    } else{
        COL <- COL.SCALE
    }
    SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), i)
    COL <- COL[SOM.hc]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight')
    add.cluster.boundaries(SOM, SOM.hc, lwd=10)
    invisible(dev.off())
}
cat("Done\n")

cat("Creating Silhouette Images...\n")
FOLDER <- paste(OUT, "/Cluster_Number/Silhouette/", sep='')
dir.create(FOLDER, showWarnings=FALSE)
#Manual Average silhouette method
SIL <- NULL
for(i in 2:50){
    png(paste(FOLDER, "/", OUT, "_Silhouette-", i, ".png", sep=''), width=2200, height=2200)
        par(cex=4, lwd=3)
        SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), i)
        sil = silhouette(SOM.hc, dist(SOM$codes[[1]]))
        plot(sil)
        abline(v=mean(sil[,3]), col='red', lty=3, lwd=2)
    invisible(dev.off())
    SIL <- c(SIL, mean(sil[,3]))
}
png(paste(FOLDER, "/", OUT, "_Silhouette-Score.png", sep=''), width=2300, height=1700)
    par(cex=4, lwd=2)
    plot(c(2:50), SIL, type='b', pch=19, lwd=2, xlab="Number of clusters", ylab='Average silhouettes', axes = FALSE)
    axis(side = 1, at = seq(0, 50, by=5), lwd=2)
    axis(side = 2, lwd=2)
    box()
invisible(dev.off())
cat("Done\n")

### POPULATION PER-REPLICA
cat("Creating Per-replica Population Images...\n")
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

SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)
for(SIM in sprintf("%04d", c(1:NREP))){
    N <- as.numeric(SIM)
    png(paste(FOLDER, "/", OUT, "_Rep_", SIM, "_population.png", sep=''), width=2200, height=2200)
    # par(cex=3.5, mfrow=c(5,4))
    par(cex=4)
    COL <- COL.SCALE[SOM.hc]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main=sprintf("Replica %04d", N))
    add.cluster.boundaries(SOM, SOM.hc, lwd=15)
    POP <- NULL
    for(NEURON in 1:nrow(SOM$grid$pts)){
        POP <- c(POP, length(which(SOM$unit.classif[STARTS[N]:ENDS[N]]==NEURON)))
    }
    POP <- POP/max(POP)
    points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=POP*MAGNIF*1.2, col="black")
    points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=POP*MAGNIF, col="white")
    invisible(dev.off())
}

#Plots with the logarithm of the population.
dir.create(paste(OUT, "/Population/Log", sep=''), showWarnings=FALSE)
SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)
for(SIM in sprintf("%04d", c(1:NREP))){
    N <- as.numeric(SIM)
    png(paste(FOLDER, "/Log/", OUT, "_Rep_", SIM, "_population.png", sep=''), width=2200, height=2200)
    # par(cex=3.5, mfrow=c(5,4))
    par(cex=4)
    COL <- COL.SCALE[SOM.hc]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main=sprintf("Replica %04d", N))
    add.cluster.boundaries(SOM, SOM.hc, lwd=15)
    POP <- NULL
    for(NEURON in 1:nrow(SOM$grid$pts)){
        POP <- c(POP, length(which(SOM$unit.classif[STARTS[N]:ENDS[N]]==NEURON)))
    }
    POP <- POP/length(SOM$unit.classif[STARTS[N]:ENDS[N]])
    POP <- log(POP)
    M <- min(POP[is.finite(POP)])
    POP <- POP-(M-1)
    POP[which(POP=="-Inf")] <- 0
    POP <- POP/max(POP)
    points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=POP*MAGNIF*1.2, col="black")
    points(SOM$grid$pts[,1], SOM$grid$pts[,2], pch=16, cex=POP*MAGNIF, col="white")
    invisible(dev.off())
}
cat("Done\n")


### PATHS for each replica
cat("Creating Per-replica Paths Images...\n")
FOLDER <- paste(OUT, "/Paths", sep='')
dir.create(FOLDER, showWarnings=FALSE)
#Stride value to depict the path
STRIDE <- 10
SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)
for(SIM in sprintf("%04d", c(1:NREP))){
    N <- as.numeric(SIM)
    png(paste(FOLDER, "/", OUT, "_Rep_", SIM, "_Path.png", sep=''), width=2000, height=2000)
    # par(cex=3.5, mfrow=c(5,4))
    par(cex=4)
    COL <- COL.SCALE[SOM.hc]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main=sprintf("Replica %04d", N))
    add.cluster.boundaries(SOM, SOM.hc, lwd=15)
    X <- NULL
    Y <- NULL
    BWR <- colorRampPalette(c("blue", "white", "red"))
    for(i in seq(STARTS[N], ENDS[N], by=STRIDE)){
        u <- SOM$unit.classif[i]
        X <- c(X, SOM$grid$pts[u,1])
        Y <- c(Y, SOM$grid$pts[u,2])
    }
    points(X,Y, pch=16, cex=30/DIM)
    points(X,Y, pch=16, col=BWR(length(X)), cex=24/DIM)
    lines(X,Y, pch=16, lwd=5)
    invisible(dev.off())
}
cat("Done\n")

### CLUSTERING THE PATHs (Only good if replicas of the same length (Steered))
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

if(PATH_CLUS=="dependent"){
    if(length(REP)==1){
        cat("Replicas of same length, Clustering Paths time-dependent...\n")
        #Function to compute the distance between two paths
        dist_path <- function(A, B, GRID){
            #A and B are two vectors of the same length containing the path through neurons while grid is the SOM$grid$pts
            D <- NULL
            for(i in 1:length(A)){
                D <- c(D, dist(rbind(GRID[A[i],], GRID[B[i],]), method="euclidean"))
            }
            return(sum(D)/length(A))
        }
    } else{
        cat("WARNING: Specified replicas of different length, using time-independent clustering...\n")
        PATH_CLUS="independent"
    }
}
if(PATH_CLUS=="independent"){
    cat("Clustering Paths time-independent...\n")
    #Function to compute the distance between two paths
    dist_path <- function(A, B, GRID){
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
}

if(NREP > 2){
   #Compute the distance matrix
   MAT <- matrix(0, ncol=NREP, nrow=NREP)
   for(i in 1:NREP){
       for(j in 1:NREP){
           MAT[i,j] <- dist_path(na.omit(steer_PATHS[,i]), na.omit(steer_PATHS[,j]), SOM$grid$pts)
       }
   }
   path.clust <- hclust(as.dist(MAT), method="complete")
   W=70*NREP
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
   cat("Done\n")
}



##################### COMPUTE TRANSITION MATRIX ################
if(GRAPH==TRUE){
    print("Generating Network analysis...")
    FOLDER <- paste(OUT, "/Transition_Matrix", sep='')
    dir.create(FOLDER, showWarnings=FALSE)
    
    #Create vector containing the classification (Neuron) assigned to each Frame
    CLASSIF <- SOM$unit.classif
    #Replace the last frame with 0 to avoid the passage from last frame to the first
    CLASSIF[STARTS] <- 0
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
    write.table(TRANS, paste(FOLDER, "/", OUT, "_Transition_Matrix.dat", sep=''))
    
    #Create transition network
    d <- NULL
    for(i in 1:dim(TRANS)[1]){
        for(j in 1:dim(TRANS)[2]){
            if(TRANS[i,j]>0){
                d <- rbind(d, c(i, j, TRANS[i,j]))
            }
        }
    }
    write.table(d, paste(FOLDER, "/", OUT, "_Transition_Network.dat", sep=''), col.names=FALSE, row.names=FALSE)
    #Remove elements on the diagonal
    D <- NULL
    for(i in 1:nrow(d)){
        if(d[i,1]!=d[i,2]){
            D <- rbind(D, d[i,])            
        }
    }
    write.table(D, paste(FOLDER, "/", OUT, "_Transition_Network-noDiagonal.dat", sep=''), col.names=FALSE, row.names=FALSE)

}
    
########################### Create NDX for the extraction of frames within neurons #########################

cat("Generating file for the Extraction of Neuron Centroids...\n")
REPRESENTATIVE <- NULL
for(i in 1:nrow(SOM$grid$pts)){
    SET <- which(SOM$unit.classif==i)
    if(length(SET)==0){
        REPRESENTATIVE <- c(REPRESENTATIVE, 0)
    } else{
        REPRESENTATIVE <- c(REPRESENTATIVE, SET[which(SOM$distances[SET]==min(SOM$distances[SET]))])
    }
}
FOLDER <- paste(OUT, "/Neurons", sep='')
dir.create(FOLDER, showWarnings=FALSE)

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

for(i in 1:max(SOM$unit.classif)){
    #Create the sub-directory of the Neuron if it does not exist
    if(dir.exists(paste(FOLDER, sprintf("/Neuron_%04d", i), sep=''))==FALSE){
        dir.create(paste(FOLDER, sprintf("/Neuron_%04d", i), sep=''))
    }
    write(sprintf("[ NEURON_%04d ]", i), file=paste(FOLDER, sprintf("/Neuron_%04d/NEURON.ndx", i), sep=''))
    NUMs <- which(SOM$unit.classif==i)
    for (j in 0:(as.integer(length(NUMs)/15)-1)){
        L <- NUMs[(j*15)+1]
        for(k in ((j*15)+2):((j+1)*15)){
            L <- paste(L, NUMs[k], sep=' ')    
        }
        write(L, file=paste(FOLDER, sprintf("/Neuron_%04d/NEURON.ndx", i), sep=''), append=TRUE)        
    }
    L <- NULL
    if((j+1)*15<length(NUMs)){
        for(k in NUMs[c((((j+1)*15)+1):length(NUMs))]){
            L <- paste(L, k, sep=' ')
        }
    }
    write(L, file=paste(FOLDER, sprintf("/Neuron_%04d/NEURON.ndx", i), sep=''), append=TRUE)
    write(paste("cd ", sprintf("Neuron_%04d/\n", i), "gmx trjconv -f $SIM -s $GRO -sub NEURON.ndx -o NEURON.xtc <<EOC\n0\nEOC\ncd ..\n",     sep=''), file=paste(FOLDER, "/Extract_Frames.sh", sep=''), append=TRUE)
}


### SOM with Neuron Numbering and cluster labels
SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)
H <- 2300+145*as.integer(NCLUS/4)
png(paste(FOLDER, "/", OUT, "_SOM-Numbering.png", sep=''), width=2200, height=H)
    par(cex=4)
    COL <- COL.SCALE[SOM.hc]
    plot(SOM, type = "mapping", bgcol=COL, col=rgb(0,0,0,0), shape='straight', main="SOM Numbering")
    add.cluster.boundaries(SOM, SOM.hc, lwd=15)
    X <- NULL
    Y <- NULL
    for(i in c(1:nrow(SOM$grid$pts))){
        X <- SOM$grid$pts[i,1]
        Y <- SOM$grid$pts[i,2]
        text(x=X, y=Y, labels=i, cex=(13/DIM), xpd=TRUE)
    }
    LET <- NULL
    for(i in LETTERS){
        for(j in LETTERS){
        LET <- c(LET, paste(i, j, sep=''))
        }
    }
    if(NCLUS < length(LETTERS)){
        LEG_LAB <- paste("Cluster ", LETTERS, sep=" ")[1:NCLUS]
        CX <- 1.24
    } else{
        LEG_LAB <- paste("Cluster ", LET, sep=" ")[1:NCLUS]
        CX <- 1.24
    }
    NR=ceiling(NCLUS/4)
    MyOrder = matrix(1:(NR*4), nrow=NR, ncol=4, byrow=T)
    
    MyBorders = rep("black", NR*4)
    MyBorders[MyOrder > NCLUS] <- NA
    legend(x=0.4, y=0, legend=LEG_LAB[MyOrder], fill=COL.SCALE[1:NCLUS][MyOrder], ncol=4, xpd=TRUE, cex=CX, bty="n", border=MyBorders) #or pch=22
invisible(dev.off())

cat("Done\n")


############################ CLUSTER CENTRO-TYPE #################################

cat("Generating file for the Extraction of Cluster Centroids...\n")
FOLDER <- paste(OUT, "/Cluster_Centroid", sep='')
dir.create(FOLDER, showWarnings=FALSE)

# write("Cluster\tFrame", file=paste(FOLDER, "/Cluster_Centroid.dat", sep=''))

SOM.hc <- cutree(hclust(dist(SOM$codes[[1]], method=DIST_CLUS), method=CLUS_METHOD), NCLUS)
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

COMPUTE_DISTANCE <- function(V1, V2){
    return(as.numeric(dist(rbind(V1, V2))))
}
if(NCLUS < length(LETTERS)){
    CLUS_NAME <- paste("Cluster", LETTERS, sep="_")[1:NCLUS]
} else{
    CLUS_NAME <- paste("Cluster", LET, sep="_")[1:NCLUS]
}
CENTROID <- sapply(unique(SOM.hc), clust.centroid, SOM$codes[[1]], SOM.hc)


Centrotypes <- NULL
invisible(suppressWarnings(file.remove(paste(FOLDER, "/Cluster_Representatives.dat", sep=''))))
invisible(file.create(paste(FOLDER, "/Cluster_Representatives.dat", sep='')))
for(i in 1:ncol(CENTROID)){
    FR <- which(SOM.hc==i)
    DistCentroid <- apply(SOM$codes[[1]], 1, COMPUTE_DISTANCE, V2=CENTROID[,i])
    Centrotypes <- which(DistCentroid==min(DistCentroid[FR]))
    write(paste(CLUS_NAME[i], ": Neuron ",  Centrotypes, sep=''), file=paste(FOLDER, "/Cluster_Representatives.dat", sep=''), append=TRUE)
#     write(Centrotypes, file=paste(FOLDER, "/Cluster_Representatives.dat", sep=''), append=TRUE)
#     write(paste(i, Centrotypes, sep="\t"), file=paste(FOLDER, "/Cluster_Centroid.dat", sep=''), append=TRUE)
}
cat("Done\n")


#Centrotypes <- NULL
#invisible(suppressWarnings(file.remove(paste(FOLDER, "/CLUSTERS.ndx", sep=''))))
#invisible(file.create(paste(FOLDER, "/CLUSTERS.ndx", sep='')))
#for(i in 1:ncol(CENTROID)){
#    FR <- which(SOM.hc==i)
#    DistCentroid <- apply(SOM$codes[[1]], 1, COMPUTE_DISTANCE, V2=CENTROID[,i])
#    Centrotypes <- which(DistCentroid==min(DistCentroid[FR]))
#    write(paste("[ ", CLUS_NAME[i], " ]", sep=''), file=paste(FOLDER, "/CLUSTERS.ndx", sep=''), append=TRUE)
#    write(Centrotypes, file=paste(FOLDER, "/CLUSTERS.ndx", sep=''), append=TRUE)
#     write(paste(i, Centrotypes, sep="\t"), file=paste(FOLDER, "/Cluster_Centroid.dat", sep=''), append=TRUE)
#}
#cat("Done\n")
#
#write("#!/bin/bash", file=paste(FOLDER, "/Extract.sh", sep=''))
#write("SIM=", file=paste(FOLDER, "/Extract.sh", sep=''), append=TRUE)
#write("GRO=\n", file=paste(FOLDER, "/Extract.sh", sep=''), append=TRUE)
#write(paste("gmx trjconv -f $SIM -s $GRO -sub CLUSTERS.ndx -o CLUSTERS.xtc <<EOC\n0\nEOC", sep=''), file=paste(FOLDER, "/Extract.sh", sep=''), append=TRUE)

save.image(paste(OUT, "/", OUT, "_SOM-Project.Rdata", sep=''))
