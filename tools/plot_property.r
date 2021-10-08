#!/usr/bin/Rscript

library(kohonen)

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
      Usage:
      R CMD BATCH plot_property.r [SOM R file] [property file]

      The property file is a single column file, with the property for each frame
      of the simulation in every row.
      If multiple replicas were analyzed by PathDetect-SOM.r, the property file is
      a single file for all the replicas concatenated.
      If a stride were applied to PathDetect-SOM.r, the same stride should be applied
      to the trajectory before computing the property value, so that the number of
      rows of the property file is equal to the number of frame parsed by
      PathDetect-SOM.r \n\n')
 
  q(save="no")
}

load(args[1])
P <- read.table(args[2])
P <- as.vector(unlist(P))

if(length(SOM$unit.classif) != length(P)){
    print(paste("ERROR: SOM input frames were ", length(SOM$unit.classif), " while rows of property file are ", length(P), sep=''))
    q(save="no")
}

AVG_NEUR_P <- NULL
for(i in 1:nrow(SOM$grid$pts)){
    AVG_NEUR_P <- c(AVG_NEUR_P, mean(P[which(SOM$unit.classif==i)]))
}

png("SOM_property_map-lab.png", width=2500, height=2500)
par(cex=4)
plot(SOM, type = "property", property=AVG_NEUR_P, shape='straight', palette.name=colorRampPalette(c("blue", "white", "red")))
for(i in 1:nrow(SOM$grid$pts)){
    if(is.na(AVG_NEUR_P[i])==FALSE){
        text(SOM$grid$pts[i,1], SOM$grid$pts[i,2], labels=round(AVG_NEUR_P[i], 1), cex=10/SOM$grid$xdim)
    }
}
dev.off()

png("SOM_property_map.png", width=2500, height=2500)
par(cex=4)
plot(SOM, type = "property", property=AVG_NEUR_P, shape='straight', palette.name=colorRampPalette(c("blue", "white", "red")))
dev.off()
