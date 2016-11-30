# Jurua analyses
# this script will be used to run the analyses of the Jurua top-hydro paper
# owner: Gabriel Massaine Moulatlet
# contact: gamamo@utu.fi


getwd()

# load the relevant packages

library(vegan)
library(ggplot2)

# this script starts with the preparation of table that gets the relative elevational of the transects in each
# geological formation, then come the graphics that relate environment and species

# load the environmental data ####

setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses/dados ambientais originais")
dir()


HAND  <- read.csv("hands_srtmnovo.csv")
topo  <- read.csv("srtm_topography.csv")
SRTM  <- read.csv("srtm_subunits.csv")
geoID <- read.csv("geoID.csv")
 
  # subset transects per geological formations
  solimoes <- geoID[which(geoID$surface == "Pebas"  ),"transect"]
  ica      <- geoID[which(geoID$surface == "Hills"  ),"transect"]
  terraces <- geoID[which(geoID$surface == "Terrace"),"transect"]
  geolist  <- list(ica, solimoes, terraces) # this list is important for the loops below


# table 1: relative elevation differences per geological formations ####
  
  statsminmax <- data.frame (
          aggregate(topo$topography, list(topo$TrNumber), min )[,1:2],
          aggregate(topo$topography, list(topo$TrNumber), max )[,-1 ],
          aggregate(topo$topography, list(topo$TrNumber), mean)[,-1 ],
          aggregate(SRTM$srtmnovo  , list(SRTM$tr      ), min )[,-1 ],
          aggregate(SRTM$srtmnovo  , list(SRTM$tr      ), max )[,-1 ],
          aggregate(SRTM$srtmnovo  , list(SRTM$tr      ), mean)[,-1 ],
          aggregate(HAND$handnew   , list(HAND$tr      ), min )[,-1 ],
          aggregate(HAND$handnew   , list(HAND$tr      ), max )[,-1 ],
          aggregate(HAND$handnew   , list(HAND$tr      ), mean)[,-1 ]
  )

  colnames(statsminmax) <- c("transect", "min_topo" , "max_topo", "mean_topo", "min_srtm",
                             "max_srtm", "mean_srtm", "min_hand", "max_hand" , "mean_hand")
  
  
  table1 <- matrix(data = NA, nrow = 3, ncol = 3) # creat the table 
    rownames(table1) <- c("Içá" , "Solimões", "Terraces")
    colnames(table1) <- c("Topo", "SRTM"    , "HAND"    )
   
    
    for(i in 1:3){
            table1[i,1] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_topo"]),2),
                                "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_topo"]),2),
                                "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_topo"]),2),")")
            table1[i,2] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_srtm"]),2),
                                 "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_srtm"]),2),
                                 "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_srtm"]),2),")")
            table1[i,3] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_hand"]),2),
                                 "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_hand"]),2),
                                 "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_hand"]),2),")")
    }
    

# graph relative topo x relative srtm and hand differences ####
  # first create a table with the differences
    
  statsminmaxdiff <- data.frame(
        statsminmax$transect,
        statsminmax$max_topo - statsminmax$min_topo,
        statsminmax$max_srtm - statsminmax$min_srtm,
        statsminmax$max_hand - statsminmax$min_hand
        )
  colnames(statsminmaxdiff) <- c("transect", "diff_topo", "diff_srtm", "diff_hand")
  
    plot(statsminmaxdiff$diff_topo, statsminmaxdiff$diff_srtm, pch=19, col="gray", xlab = "relative topographic differences",
       ylab = "relative remote sensing differences", ylim = c(0,40), xlim = c(0,40))
          r <- lm(statsminmaxdiff$diff_srtm~statsminmaxdiff$diff_topo)
          predict.lm(r)
          lines(x=c(1,40), y=c(min(predict.lm(r)),max(predict.lm(r))),col="gray")
    par(new=T)
    plot(statsminmaxdiff$diff_topo, statsminmaxdiff$diff_hand, pch=19, col="red", ylim = c(0,40), xlim = c(0,40), xlab="",
       ylab = "")
          abline(lm(statsminmaxdiff$diff_topo~statsminmaxdiff$diff_hand), col="red")
    legend("bottomright", legend = c("SRTM", "HAND"), text.col = c("gray","red"),bty = "n", cex=2)
    
    
  
  

