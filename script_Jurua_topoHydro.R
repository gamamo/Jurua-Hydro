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
  colnames(HAND)[2:3] <- c("TrNumber", "subunit") #change names of HAND columns
  HAND$subunit <- (HAND$subunit+12.5)/25          #change values of subunit column. They are in a different format
  
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
          aggregate(HAND$handnew   , list(HAND$TrNumber), min )[,-1 ],
          aggregate(HAND$handnew   , list(HAND$TrNumber), max )[,-1 ],
          aggregate(HAND$handnew   , list(HAND$TrNumber), mean)[,-1 ]
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
    print(table1)
    
    getwd()       # sabe in the folder "outputs
    setwd("..")   # one folder up
    
    file <- paste(getwd(),"outputs","table1.csv",sep="/")
    write.csv(table1, file)
    
# graph 1: relative topo x relative srtm and hand differences ####
  # first create a table with the differences
  # then create the graphic and then export as pdf
    
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

# importing species data, run ordinations and plot against environmental variables ####
  
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
        names(species25list) <- c("ferns", "zings", "palms")
        
        
       # call the environmental data
       # each plant object has a different number of rows. It has to be adequated in the environmental data: use the join::plyr 
        
        HAND
        SRTM
        topo
        geoID
        geolist
        
        
       # calculate species optima and toleraces to the HAND gradient
        library(rioja)
        library(palaeoSig)
        
        listrelabutotal     <- list()
        listrelabutotalHAND <- list()
        
        for(i in seq(length(species25list))){
          sprabd <- decostand(species25list[[i]][-c(1:2)], method = "total",1)      #to calculate relative abundances
          listrelabutotal[[i]]     <- sprabd     #save the output
          listrelabutotalHAND[[i]] <- join(species25list[[i]][c(1:2)], HAND, by = c("TrNumber","subunit"), type = "left", match = "all") # and merge with envi
        }
          
        # graphic 2: species tolerances along the HAND gradient ####
        
        graphnamesopt <- names(species25list)
        par(mfrow = c(1,3))
        for(i in seq(length(listrelabutotal))){
          fit_25 <- WA(listrelabutotal[[i]], log10(1+listrelabutotalHAND[[i]][,"handnew"]), tolDW=TRUE) 
          coef_25 <- as.data.frame(coefficients(fit_25))
          par(mar=c(5.1,6.1,2.1,1.1))
          centipede.plot(fit_25, main=graphnamesopt[[i]])

        }


       # run the NMDS ordinations
       # first calculate the relative abundances
       # then run the NMDS and store the first 2 axis of each plant group in a separate list
        library(vegan)
        library(plyr)
    
        
        listnmdstotal <- list()
        
        for (i in seq(length(species25list))){
          
          dist.ab <- vegdist(decostand(species25list[[i]][-c(1:2)], method = "total",1), method = "bray")
          mds.ab  <- monoMDS(dist.ab, y = cmdscale(dist.ab, k=2),k = 2, model = "global", threshold = 0.8, maxit = 200, 
                             weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
          listnmdstotal[[i]] <- cbind(species25list[[i]][c(1:2)], scores(mds.ab)[,c(1:2)])
          
        }
        
      
      #  graphic 3: total NMDS1 x HAND ####
        
        HANDsubsettotal <- list()
        
        for (i in seq(length(listnmdstotal))){
          HANDsubsettotal[[i]] <- join(listnmdstotal[[i]], HAND, by = c("TrNumber","subunit"), type = "left", match = "all")
        }
    
        
        graphnames <- c("Ferns","Zings","Palms")
        par(mfrow = c(1,3))
        for(i in seq((listnmdstotal))){
          plot(HANDsubsettotal[[i]][,"handnew"],listnmdstotal[[i]][,"MDS1"], xlab = "HAND", ylab = "NMDS",
               pch=19, col="gray",main = graphnames[i],ylim = c(-3,3))
        }
    
    
      # graphic 4: NMDS1 x HAND per geoformation ####
      # NMDS for each geoformation has to be generated again  
      # they will be stored in a list
        
        library(vegan)
        
        listnmdsgeo <- list()
        
        for(g in seq(length(geolist))){                 #loop over the geo formations
              for (i in seq(length(species25list))){    #loop over the plant groups
    
          dist.ab <- vegdist(decostand(species25list[[i]][which(match(species25list[[i]][,1],geolist[[g]])>0),-c(1:2)]   #run nmds
                                       , method = "total",1), method = "bray")
          mds.ab  <- monoMDS(dist.ab, y = cmdscale(dist.ab, k=2),k = 2, model = "global", threshold = 0.8, maxit = 200,  #get scores
                             weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
          trandsubinfo <- species25list[[i]][which(match(species25list[[i]][,1],geolist[[g]])>0),c(1:2)]                 #save tr and subunit info
          
         if (g==1)
         {listnmdsgeo[[i]]  <- cbind(trandsubinfo,scores(mds.ab)[,c(1:2)])}   # store the results in the right list slot
         if (g==2)
         {listnmdsgeo[[i+3]]<- cbind(trandsubinfo,scores(mds.ab)[,c(1:2)])}   # store the results in the right list slot
         if (g==3)  
         {listnmdsgeo[[i+6]]<- cbind(trandsubinfo,scores(mds.ab)[,c(1:2)])}   # store the results in the right list slot
          }
        }
        names(listnmdsgeo) <- c("fern_ica","zing_ica","palm_ica",             # give names for the list elements
                                "fern_soli","zing_soli","palm_soli",
                                "fern_ter","zing_ter","palm_ter")
        
        
        
        # NMDS1 x per subunit per geo formation. Merge the listnmdsgeo with HAND to get the exaclty rows for each species group
        # in each geoformation
        library(plyr)
    
        HANDsubsetgeo <- list()
        
        for (i in seq(length(listnmdsgeo))){
          HANDsubsetgeo[[i]] <- join(listnmdsgeo[[i]], HAND, by = c("TrNumber","subunit"), type = "left", match = "all")
        }
        
        names(HANDsubsetgeo) <- c("fern_ica","zing_ica","palm_ica",             # give names for the list elements
                                  "fern_soli","zing_soli","palm_soli",
                                  "fern_ter","zing_ter","palm_ter")
        
        # graphic 5: NMDS1 x per subunit per geo formation ####
        #obs: there is weird nmds1 value for ferns_terrace: check!!!
        
        graphnames2 <- names(listnmdsgeo)
        par(mfrow = c(3,3))
        
        for (i in seq(length(listnmdsgeo))) {
          
          plot(HANDsubsetgeo[[i]][,"handnew"],listnmdsgeo[[i]][,"MDS1"],xlab = "HAND", ylab = "NMDS",
               pch = 19, col = "gray", main = graphnames2[i])
        }
        
    
# import the moisture indexes and plot them against environmental variables ####
        
        getwd()
        setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses")
        
        HAND
        SRTM
        topo
        geoID
        geolist
        
        fernsMI <- read.csv("MIferns_persubunit.csv")
        zingsMI <- read.csv("MIzing_persubunit.csv")
        palmsMI <- read.csv("MIpalm_persubunit.csv")
        
        # creat a list with all the MIs
        
        MIlist <- list(fernsMI, zingsMI, palmsMI)
        names(MIlist) <- c("fernsMI","zingsMI", "palmsMI")
        
        # graphic 6: plot MI x HAND total ####
        library(plyr)
        
        HANDsubsetMI <- list()
        
        for (i in seq(length(MIlist))){
          HANDsubsetMI[[i]] <- join(MIlist[[i]], HAND, by = c("TrNumber","subunit"), type = "left", match = "all")
        }
        
        
        graphnames <- c("Ferns","Zings","Palms")
        par(mfrow = c(1,3))
        for(i in seq((MIlist))){
          plot(HANDsubsetMI[[i]][,"handnew"],MIlist[[i]][,"MI"], xlab = "HAND", ylab = "MI",
               pch=19, col="gray",main = graphnames[i],ylim = c(0,2))
        }
        
        # plot MI x HAND per geo formation
        # first creat a list of MI per species per geo
        
        MIlistgeo <- list()
    
        for(g in seq(geolist)){
           for(i in seq(MIlist)){
            subsetMIlist <- MIlist[[i]][which(match(MIlist[[i]][,"TrNumber"],geolist[[g]])>0),]
            if (g==1)
            {MIlistgeo[[i]]  <- subsetMIlist}
            if (g==2)
            {MIlistgeo[[i+3]]<- subsetMIlist}
            if (g==3)  
            {MIlistgeo[[i+6]]<- subsetMIlist}    
          }
        }
        names(MIlistgeo) <- c("MIferns_ica", "MIzings_ica", "MIpalms_ica",
                                "MIferns_sol", "MIzings_sol", "Mipalms_sol",
                                "MIferns_ter", "MIzings_ter", "MIpalms_ter")
    
        
        # graphic 7: MI x HAND per geo formation 
        # first merge HAND values with the MIlistgeo
    
        library(plyr)
        
        HANDsubsetgeoMI <- list()
        
        for (i in seq(length(MIlistgeo))){
          HANDsubsetgeoMI[[i]] <- join(MIlistgeo[[i]], HAND, by = c("TrNumber","subunit"), type = "left", match = "all")
        }
        
        names(HANDsubsetgeoMI) <- c("fern_ica","zing_ica","palm_ica",             # give names for the list elements
                                    "fern_soli","zing_soli","palm_soli",
                                    "fern_ter","zing_ter","palm_ter")
        
        
        # graphic 8: MI x HAND per geo formation ####
        
        graphnames3 <- names(MIlistgeo)
        par(mfrow = c(3,3))
        
        for (i in seq(length(MIlistgeo))) {
          plot(HANDsubsetgeoMI[[i]][,"handnew"], MIlistgeo[[i]][,"MI"],xlab = "HAND", ylab = "MI",
               pch = 19, col = "gray", main = graphnames3[i])
        }

        
        
        ### testes
        
    HANDsubsetgeoMI[[7]][101,]
          y<- MIlistgeo[[7]][,"MI"]
          x <- HANDsubsetgeoMI[[7]][,"handnew"]
          
          cor(HANDsubsetgeoMI[[7]][,"handnew"],MIlistgeo[[7]][,"MI"])
        
        a <- data.frame(HANDsubsetgeoMI[[7]][,"handnew"], MIlistgeo[[7]][,"MI"])
        colnames(a) <- c("hand","mi")
    
    library(ggplot2)
            ggplot(a, aes(a$hand,a$mi)) + geom_point() + geom_smooth()
dev.off()
