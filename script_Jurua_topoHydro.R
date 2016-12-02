# Jurua analyses
# this script will be used to run the analyses of the Jurua top-hydro paper
# owner: Gabriel Massaine Moulatlet
# contact: gamamo@utu.fi


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
    
    getwd()
    
    file <- paste(getwd(),"outputs","table1.csv",sep="/")
    write.csv(table1, file)
    
# graph relative topo x relative srtm and hand differences ####
  # first create a table with the differences
  # then create the graphic and export to pdf  
    
  statsminmaxdiff <- data.frame(
        statsminmax$transect,
        statsminmax$max_topo - statsminmax$min_topo,
        statsminmax$max_srtm - statsminmax$min_srtm,
        statsminmax$max_hand - statsminmax$min_hand
        )
  colnames(statsminmaxdiff) <- c("transect", "diff_topo", "diff_srtm", "diff_hand")
  
#here is the command for the graphic
  
  #save in the folder "outputs"
  
  getwd()
  setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses")
  
  #run the plot
  
  file <- paste(getwd(),"outputs","figure_relative_diff.pdf",sep="/") #always check if the file is getting saved 
  mfrow <- c(1,1)
  pdf(file=file, height = 7*mfrow[1], width = 8*mfrow[2])
  par(mfrow=mfrow, mar=c(0.2,0.8,2,0.5), oma=c(2,1,.5,0.5), mgp=c(1.7,0.6,0))
  
    plot(statsminmaxdiff$diff_topo, statsminmaxdiff$diff_srtm, pch=19, col="black", xlab = "relative topographic differences",
       ylab = "relative remote sensing differences", ylim = c(0,40), xlim = c(0,40),cex.lab=1.5)
          r1 <- lm(statsminmaxdiff$diff_srtm~statsminmaxdiff$diff_topo)
          predict.lm(r1)
          lines(x=c(1,40), y=c(min(predict.lm(r1)),max(predict.lm(r1))),col="black")
    par(new=T)
    plot(statsminmaxdiff$diff_topo, statsminmaxdiff$diff_hand, pch=19, col="gray", ylim = c(0,40), xlim = c(0,40), xlab="",
       ylab = "")
          r2 <- lm(statsminmaxdiff$diff_hand~statsminmaxdiff$diff_topo)
          predict.lm(r2)
          lines(x=c(1,40), y=c(min(predict.lm(r2)),max(predict.lm(r2))),col="gray")
    legend("bottomright", legend = c(paste("SRTM","-","R2",round(summary(r1)$r.squared,2)), 
                                     paste("HAND","-","R2",round(summary(r2)$r.squared,2))),
           text.col = c("black","gray"),bty = "n", cex=1.3)
  dev.off()   

# PART 2: importing species data, run ordinations and moisture indexes ####
  
  # import species data
  # these species tables were generated in a script called "script_preparacao_sp_analises.R"
  # the files are in the folder "dados floristicos originais"
  # navigate to that directory
  
  getwd()
  setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses")

  fern25 <- read.csv("ferns25_widetableNE.csv"     , stringsAsFactors = FALSE)
        fern25 <- fern25[-which(rowSums(fern25[,-c(1:2)])=="0"),]  # delete subunits with zero occurrences
        fern25 <- fern25[-c(795,1210),] # this 2 subunits were deleted because were too weird. Have to check it again
      
  zing25 <- read.csv("zingdata_wide_gmm_v1.csv"    , stringsAsFactors = FALSE)
        #zing25 <- zing25[-which(rowSums(zing25[,-c(1:2)])=="0"),]  # delete subunits with zero occurrences      
  
  palm25 <- read.csv("palms_jurua_subunit_gmm1.csv", stringsAsFactors = FALSE) # no need to delete rows for the palms: checked before
            rowSums(palm25[,-c(1:2)])
            
    species25list <- list(fern25, zing25, palm25) # creat a list if the species data


  # run the NMDS ordinations
  # first calculate the relative abundances
  # then run the NMDS and store the first 2 axis of each plant group in a separate list
    library(vegan)
    
    listnmds <- list()
    
    for (i in seq(length(species25list))){
      
      dist.ab <- vegdist(decostand(species25list[[i]][-c(1:2)], method = "total",1), method = "bray")
      mds.ab  <- monoMDS(dist.ab, y = cmdscale(dist.ab, k=2),k = 2, model = "global", threshold = 0.8, maxit = 200, 
                         weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
      listnmds[[i]] <- scores(mds.ab)[,c(1:2)]  
      
    }

  # call the environmental data
  # each plant object has a differene number of rows. It has to be adequated in the environmental data  
    
    HAND
    SRTM
    topo
    geoID
    geolist
  
  #  graphic 1: total NMDS x HAND
    
    graphnames <- c("Ferns","Zings","Palms")
    par(mfrow = c(1,3))
    for(i in seq((listnmds))){
      plot(HAND[rownames(species25list[[i]]),"handnew"],listnmds[[i]][,1], xlab = "HAND", ylab = "NMDS",
           pch=19, col="gray",main = graphnames[i],ylim = c(-3,3))
    }

    
    
  # graphic 2: NMDS x HAND per subunit
  # NMDS for each geoformation has to be generated again  
    library(vegan)
    
    listnmdsgeo <- list()
    
    for(g in seq(length(geolist))){
          for (i in seq(length(species25list))){

      dist.ab <- vegdist(decostand(species25list[[i]][which(match(species25list[[i]][,1],geolist[[g]])>0),-c(1:2)]
                                   , method = "total",1), method = "bray")
      mds.ab  <- monoMDS(dist.ab, y = cmdscale(dist.ab, k=2),k = 2, model = "global", threshold = 0.8, maxit = 200, 
                         weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
      listnmdsgeo[[i]] <- scores(mds.ab)[,c(1:2)]  
          
      }
    }
    
    species25list[[1]][which(match(species25list[[1]][,1],geolist[[1]])>0),-c(1:2)]
    