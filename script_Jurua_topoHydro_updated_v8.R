      # Jurua analyses
      # this script will be used to run the analyses of the Jurua top-hydro paper
      # owner: Gabriel Massaine Moulatlet
      # contact: gamamo@utu.fi
      
      
      # load the relevant packages
      
      library(vegan)
      library(ggplot2)
      library(lme4)
      library(plyr)
      library(dplyr)
      library(nlme)
      library(effects)
      library(car)
      library(pROC)
      library(mvpart)
      library(rioja)
      library(palaeoSig)
      library(rpart)
      library(visreg)
      library(gamlss)
      library(raster)
      library(sp)
      library(rgdal)
      library(sjPlot)
      library(coefplot)
      library(VennDiagram)
      library(reshape2)
      library(venneuler)
      library(gridExtra)
      library(pbkrtest)
      library(mgcv)
      library(voxel)
      library(ggpubr)
      library(ggrepel)
      library(glmmTMB)
      library(MVPARTwrap)
      library(rpart.plot)

      
      # this script starts with the preparation of table that gets the relative elevational of the transects in each
      # geological formation, then come the graphics that relate environment and species
      
###################################### load the environmental data ####
      
      rm(list = ls())
      
      getwd()
      setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses/dados ambientais originais")
      #setwd("F:/workspace gabriel/Hidrologia do jurua/Analyses/dados ambientais originais") 
      dir()
      
      
      
      #load GEO list
      geoID <- read.csv("geoID.csv", stringsAsFactors = FALSE)
        geoID[which(geoID$surface=="Terrace"),"surface"] <- "River terraces"
        geoID[which(geoID$surface=="Pebas"),"surface"] <- "Solimões formation"
        geoID[which(geoID$surface=="Hills"),"surface"] <- "Içá formation"
        colnames(geoID)[1:2] <- c("trname","tr")
      
      #load solos    
      solo  <- read.csv("soil_results_GMM_v3_forR.csv",stringsAsFactors = F)
        colnames(solo)[2:3] <- c("TrNumber", "subunit")

                #solo, from 1:100 to 1:20
                newsubunitnamesolo<- matrix(seq(1,100), ncol = 20) # create a matrix with numbers to be related
                
                #run the replace looping
                for(i in seq(1,20)){
                  solo$subunit[which(solo$subunit %in% newsubunitnamesolo[,i])] <- i
                }  
      
      #load the topography data
        topography <- read.csv("topography.csv", header = TRUE, sep = ",", dec = ".", na.strings = NA)
        head(topography)
            
            #identify which transects do not have topo measurments
            topography[which(topography$sub==0 & topography$topo==0),]
            # transects 780,794,806,807,808 need to be discarted
            ex <- c(780,794,806,807,808)
            
            topography <- topography[!topography$tr %in% ex,]
        
                dist <- c(1:20 * 25)
                
                #interpolate the topography to each 25 m
                listtopointerpolation <- list()
                for (k in unique(topography$tr)){
                  topo_sel <- subset(topography, topography$tr==k)
                  topoapprox <- approx(topo_sel$sub,topo_sel$topo,dist)
                  
                  matrixtopo <- matrix(NA, ncol=3, nrow = length(topoapprox$y))
                  matrixtopo[,1] <- topoapprox$y
                  matrixtopo[,2] <- rep(k, times=length(topoapprox$y))
                  matrixtopo[,3] <- rep(1:20)
                  numerator <- which(unique(topography$tr)==k)
                  listtopointerpolation[[numerator]] <- matrixtopo
                }
                
                topo_sel <- ldply(listtopointerpolation,data.frame)
                
                colnames(topo_sel) <- c("topo","trnumber","sub")
                
      #setwd("C:/workspace gabriel/hidrologia do jurua/camilo R script june 2017")
      #setwd("F:/workspace gabriel/Hidrologia do jurua/camilo R script june 2017")
      setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/camilo R script june 2017")
        dir()
        SRTM <-raster("SRTM_JURUA_1arc_rec.tif")
        
        #Load the bilinear extractions made in R by Camilo
      hand53_ream_m0 <-  read.csv("transect_adjust_m0_ream_hand.csv",sep = ";")  
          head(hand53_ream_m0)
          hand53_ream_m0$sub <- rep(seq(1:20),length(unique(hand53_ream_m0$tr)))
      hand53_ream_01<-  read.csv("transect_adjust_01_ream_hand.csv",sep = ";")  
          head( hand53_ream_01)
        
          #coordinates(hand53_ream_m0) <- c("lat","lon")
          srtm_alt <- extract(SRTM,cbind(hand53_ream_m0$lon,hand53_ream_m0$lat),method="bilinear" )
          hand53_ream_m0 <- cbind(hand53_ream_m0,srtm_alt)
          rm(SRTM)
          
################################# prepare envi table ####
          
      tr <- seq(738,808)
          envi <- list()
          for(i in unique(tr)){
            envi[[i]]<- rep(i, times=20)
          }
          envi <- ldply(envi,data.frame)
          colnames(envi) <- "tr"

          envi$sub <- NA
          envi$LOI <- NA
          envi$pH  <- NA
          envi$Ca  <- NA
          envi$K   <- NA
          envi$Mg  <- NA
          envi$Na  <- NA
          envi$sum_of_basis <- NA
          envi$Al  <- NA
          envi$P   <- NA
          envi$surface <- NA
          envi$topo <- NA
          envi$lon <- NA
          envi$lat  <- NA
          envi$hand <- NA
          envi$area_c <- NA  
          envi$decli <- NA
          envi$drain <- NA
          envi$max_drain <- NA
          envi$soloslinearint <- NA
          envi$solostopoint <- NA
          envi$srtm_alt    <- NA
      
      #assign sub information
          for(i in unique(envi$tr)){
            envi[envi$tr==i,"sub"] <- seq(1:20)
          }
          
         
          
      #associate soil properties to envi
          for(i in unique(envi$tr)){
            for(g in unique(envi$sub)){
                if(nrow(solo[solo$TrNumber==i & solo$subunit==g,])>0){
                  envi[envi$tr==i & envi$sub==g,
                       c("LOI","pH","Ca","K","Mg","Na","sum_of_basis","Al","P")] <- solo[solo$TrNumber==i & solo$subunit==g,
                                                                                         c("LOI","pH","Ca","K","Mg","Na","sum_of_basis","Al","P")]
                }else{
                  envi[envi$tr==i & envi$sub==g,c("LOI","pH","Ca","K","Mg","Na","sum_of_basis","Al","P")] <- NA
                  }
                }
          }
          
      # interpolated sum of basis
          listsolointerpolation <- list()
          for (i in unique(solo$TrNumber)){
            solo_sel <- subset(solo, solo$TrNumber==i)
            soloapprox <- approx(solo_sel$subunit*25,solo_sel$sum_of_basis,dist,method = "linear",rule=2)
            soloapprox2 <- spline(solo_sel$subunit*25,solo_sel$sum_of_basis,xout=topo_sel[topo_sel$trnumber==k,"topo"],method = "natural")
          
            matrixsolo <- matrix(NA, ncol=4, nrow = length(soloapprox$y))
            matrixsolo[,1] <- soloapprox$y
            matrixsolo[,2] <- soloapprox2$y
            matrixsolo[,3] <- rep(i, times=length(soloapprox$y))
            matrixsolo[,4] <- rep(1:20)
            numerator <- which(unique(solo$TrNumber)==i)
            listsolointerpolation[[numerator]] <- matrixsolo
          }
          unlistsolos <- ldply(listsolointerpolation,data.frame)  
          colnames(unlistsolos) <- c("soloslinearint","solostopoint","trnumber","sub")
          
          
          for (i in unique(unlistsolos$trnumber)){
            envi[envi$tr==i, "soloslinearint"] <- unlistsolos[unlistsolos$trnumber==i,"soloslinearint"]
            envi[envi$tr==i, "solostopoint"] <- unlistsolos[unlistsolos$trnumber==i,"solostopoint"]
          }

      #associate geo surface info to envi
          for (i in unique(envi$tr)){
            envi[envi$tr==i, "surface"] <- geoID[geoID$tr==i,"surface"]
          }
          
      #associate topo to envi
          for (i in unique(topo_sel$trnumber)){
            envi[envi$tr==i, "topo"] <- topo_sel[topo_sel$trnumber==i,"topo"]
          }
          
      #associate hand data
          for (i in unique(hand53_ream_m0$tr)){
            envi[envi$tr==i, c("lon","lat","hand","area_c","decli", "drain", "max_drain","srtm_alt")] <- hand53_ream_m0[hand53_ream_m0$tr==i,
                                                                                                              c("lon","lat","hand",
                                                                                                                "area_c","decli", "drain", "max_drain",
                                                                                                                "srtm_alt")]
          }
          

      #associate Truni info   
          envi$Truni <- as.numeric(paste(envi$tr,envi$sub,sep = "")) 
      
      #colnames change
          colnames(envi)[1:2] <- c("TrNumber","subunit")
          
      #transects to exclude ####
          toexclude <- as.factor(c(743,742,745,764,765,770,775,784,785,805))
          trvector <- as.factor(envi$TrNumber)
          trvector <- trvector[!trvector %in% toexclude]
          trvector  <- as.numeric(as.vector(unique(trvector)))
          
          envi <- envi[envi$TrNumber %in% trvector,]

          
# end of envi table preparition  
                
################################################# table 1: relative elevation differences per geological formations ####
          
            # subset transects per geological formations
            solimoes <- envi[which(envi$surface == "Solimões formation"  ),"TrNumber"]
            ica      <- envi[which(envi$surface == "Içá formation"       ),"TrNumber"]
            terraces <- envi[which(envi$surface == "River terraces"     ),"TrNumber"]
            geolist  <- list(ica, solimoes, terraces) # this list is important for the loops below
                
                statsminmax <- data.frame (
                        aggregate(envi$topo, list(envi$TrNumber), min         )[,1:2],
                        aggregate(envi$topo, list(envi$TrNumber), max         )[,-1 ],
                        aggregate(envi$topo, list(envi$TrNumber), mean        )[,-1 ],
                        aggregate(envi$hand, list(envi$TrNumber), min         )[,-1 ],
                        aggregate(envi$hand, list(envi$TrNumber), max         )[,-1 ],
                        aggregate(envi$hand, list(envi$TrNumber), mean        )[,-1 ],
                        aggregate(envi$sum_of_basis, list(envi$TrNumber), min ,na.rm=TRUE )[,-1 ],
                        aggregate(envi$sum_of_basis, list(envi$TrNumber), max ,na.rm=TRUE )[,-1 ],
                        aggregate(envi$sum_of_basis, list(envi$TrNumber), mean,na.rm=TRUE)[,-1 ],
                        aggregate(envi$decli, list(envi$TrNumber), min        )[,-1 ],
                        aggregate(envi$decli, list(envi$TrNumber), max        )[,-1 ],
                        aggregate(envi$decli, list(envi$TrNumber), mean       )[,-1 ],
                        aggregate(envi$srtm_alt, list(envi$TrNumber), min        )[,-1 ],
                        aggregate(envi$srtm_alt, list(envi$TrNumber), max        )[,-1 ],
                        aggregate(envi$srtm_alt, list(envi$TrNumber), mean       )[,-1 ],
                        aggregate(envi$soloslinearint, list(envi$TrNumber), min        )[,-1 ],
                        aggregate(envi$soloslinearint, list(envi$TrNumber), max        )[,-1 ],
                        aggregate(envi$soloslinearint, list(envi$TrNumber), mean       )[,-1 ]
                )
              
                colnames(statsminmax) <- c("transect", "min_topo" , "max_topo", "mean_topo", "min_hand", "max_hand" , "mean_hand",
                                           "min_cation", "max_cation" , "mean_cation","min_slope", "max_slope" , "mean_slope",
                                           "min_srtm", "max_srtm" , "mean_srtm","min_soloint", "max_soloint" , "mean_soloint")
                
                
                table1 <- matrix(data = NA, nrow = 3, ncol = 6) # creat the table 
                  rownames(table1) <- c("Içá formation" , "Solimões formation", "River terraces")
                  colnames(table1) <- c("Field Topography" , "Soil Moisture (HAND)", "Soil Fertility (Cation concentration)", "Slope",
                                        "SRTM","Cation concentration (interpolated values)")
                 
                  
                  for(i in 1:3){
                          table1[i,1] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_topo"]  ,na.rm = TRUE),2),
                                              "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_topo"],na.rm = TRUE),2),
                                              "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_topo"],na.rm = TRUE),2),")")
                          
                          table1[i,2] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_hand"]   ,na.rm = TRUE),2),
                                               "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_hand"],na.rm = TRUE),2),
                                               "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_hand"],na.rm = TRUE),2),")")
                          
                          table1[i,3] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_cation"]   ,na.rm = TRUE),2),
                                               "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_cation"],na.rm = TRUE),2),
                                               "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_cation"],na.rm = TRUE),2),")")
                          
                          table1[i,4] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_slope"]   ,na.rm = TRUE),2),
                                               "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_slope"],na.rm = TRUE),2),
                                               "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_slope"],na.rm = TRUE),2),")")
                          
                          table1[i,5] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_srtm"]   ,na.rm = TRUE),2),
                                               "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_srtm"],na.rm = TRUE),2),
                                               "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_srtm"],na.rm = TRUE),2),")")
                          
                          table1[i,6] <- paste(round(mean(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"mean_soloint"]   ,na.rm = TRUE),2),
                                               "(", round(min(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"min_soloint"],na.rm = TRUE),2),
                                               "-", round(max(statsminmax[which(match(statsminmax$transect,geolist[[i]])>0),"max_soloint"],na.rm = TRUE),2),")")
                  }
                  print(table1)
                  
                  setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses")
                  #setwd("F:/workspace gabriel/Hidrologia do jurua/Analyses") 
                  
                  write.csv(table1, "RESUtable2.csv",row.names = TRUE)
                  
#end of table 1
                  
############################################################  calculate differences wihin transects ####      

        # first create a table with the differences
        # then create a dicionary that contains the order of topogaphic differences  
          
              statsminmaxdiff <- data.frame(
                    statsminmax$transect,
                    statsminmax$max_topo - statsminmax$min_topo,
                    statsminmax$max_hand - statsminmax$min_hand
                    )
              colnames(statsminmaxdiff) <- c("transect", "diff_topo", "diff_hand")
              statsminmaxdifford <- statsminmaxdiff[order(statsminmaxdiff$diff_topo,decreasing = TRUE),]
              
              dici_topodiff <- data.frame(statsminmaxdifford$transect,statsminmaxdifford$diff_topo)
                    colnames(dici_topodiff) <- c("TrNumber","topo")
                    dici_topodiff$topodifford <- as.numeric(seq(1:length(dici_topodiff[,1])))
              
# end of table diff within transects              

########################################################## importing species data and prepare a list with envi of each species group ####
        
            # import species data
            # these species tables were generated in a script called "script_preparacao_sp_analises.R"
            # the original files are in the folder "dados floristicos originais"
            
            #setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses")
          
            fern25 <- read.csv("ferns25_widetableWE.csv"     , stringsAsFactors = FALSE)
                  fern25 <- fern25[-which(rowSums(fern25[,-c(1:2)])=="0"),]  # delete subunits with zero occurrences
                  fern25 <- fern25[fern25$TrNumber %in% trvector,]
                  
                  
                
            zing25 <- read.csv("zingdata_wide_gmm_v1.csv"    , stringsAsFactors = FALSE)
                  zing25 <- zing25[zing25$TrNumber %in% trvector,]
                  colnames(zing25)
            

            palm25 <- read.csv("palms_jurua_subunit_gmm1.csv", stringsAsFactors = FALSE) # no need to delete rows for the palms: checked before
                  palm25 <- palm25[palm25$TrNumber %in% trvector,]
                  
            melas25 <- read.csv("Mel_Jurua_5x25m_table_gmm.csv", stringsAsFactors = FALSE)    
                  melas25 <- melas25[melas25$TrNumber %in% trvector,]
                    temp <- melas25[,-c(1:2)]
                    which(rowSums(temp)!=0)
                  melas25 <- melas25[which(rowSums(temp)!=0),]
                 
                  colSums(melas25) 
                 
                  
            # creat a list of the species data to GLMMS             
              species25list <- list(fern25, zing25, palm25,melas25) 
              names(species25list) <- c("Ferns", "Zingiberales", "Arecaceae","Melastomataceae")
              
              
              
            #how many species?
              length(fern25[,-c(1:2)])
              length(zing25[,-c(1:2)])
              length(palm25[,-c(1:2)])
              length(melas25[,-c(1:2)])

            # join all species into a single dataframe
              
              join1 <- join(fern25, zing25 ,by = c("TrNumber","subunit"), type = "inner", match = "all")
              join2 <- join(join1 , palm25 ,by = c("TrNumber","subunit"), type = "inner", match = "all")
              join3 <- join(join2 , melas25,by = c("TrNumber","subunit"), type = "inner", match = "all")

              speciesall <- join3
              speciesall_tomodel <- speciesall[,-c(1:2)]
              speciesall_tomodel_pa  <- decostand(speciesall_tomodel, method = "pa",   1)   #to transform to PA
              
            # join species and envi in a dataframe
              
              speciesenvi <- join(speciesall, envi,by = c("TrNumber","subunit"), type = "inner", match = "all" )
              
              dd <- cbind(speciesenvi$Adian.term,speciesenvi$surface)
              dd[order(dd[,1],decreasing = T),]
              
              #select some ferns and zings as a shapefile to send to camilo
              
              fernstocamilo <- decostand(fern25[,c("Adian.peti","Lomag.guia")],method="pa",1)
              fernstocamilo <- cbind(fern25[,c("TrNumber","subunit")],fernstocamilo)
              fernstocamilo <- join(fernstocamilo,envi,by = c("TrNumber","subunit"), type = "left", match = "all")
                  fernstocamilo <- fernstocamilo[,c("TrNumber","subunit","Adian.peti","Lomag.guia","lon","lat")]
                  fernstocamilo <- fernstocamilo[which(fernstocamilo$lon!="NA"),]
                  coordinates(fernstocamilo) <- ~lon+lat
                  writeOGR(fernstocamilo,".","ferns_camilo2",driver ="ESRI Shapefile" )
                
              zingscamilo <- decostand(zing25[,c("cala.zing","cala.mica")],method="pa",1)
              zingscamilo <- cbind(zing25[,c("TrNumber","subunit")],zingscamilo)
              zingscamilo <- join(zingscamilo,envi,by = c("TrNumber","subunit"), type = "left", match = "all")    
                  zingscamilo <- zingscamilo[,c("TrNumber","subunit","cala.zing","cala.mica","lon","lat")]
                  zingscamilo <-  zingscamilo[which( zingscamilo$lon!="NA"),]
                  coordinates( zingscamilo) <- ~lon+lat
                  writeOGR( zingscamilo,".","zings_camilo2",driver ="ESRI Shapefile" )
              
              
                  
            # call the environmental data
            # each plant object has a different number of rows. It has to be adequated in the environmental data: use the join::plyr 
            # calculate species optima and toleraces to the HAND gradient
              
              listrelaabutotal <- list() #only abu
              listpatotal      <- list() #only pa
              listspENVI       <- list() #join abu and envi
              listspENVIpa     <- list() #join pa  and envi
              
              
              for(i in seq(length(species25list))){
                sprela <- decostand(species25list[[i]][-c(1:2)], method = "total",1)   #to calculate relative abundances
                sppa   <- decostand(species25list[[i]][-c(1:2)], method = "pa",   1)   #to transform to PA
                
                listrelaabutotal[[i]] <- sprela
                listpatotal[[i]]      <- vegdist(sppa,method="bray")
                listspENVI[[i]]       <- join(species25list[[i]], envi, by = c("TrNumber","subunit"), type = "left", match = "all") # and merge with envi
                listspENVIpa[[i]]     <- join(cbind(species25list[[i]][c(1:2)],sppa), envi, by = c("TrNumber","subunit"), type = "left", match = "all") # and merge with envi
              }

#end of species importing and envi matching


###################################################################### Run ordinations e prepare the table for the analysis #####
      
      
             # run the NMDS ordinations
             # first calculate the relative abundances
             # then run the NMDS and store the first 2 axis of each plant group in a separate list
       
              listnmdsab <- list()
              listnmdspa <- list()
              
              for (i in seq(length(species25list))){
                dist.ab <- vegdist(decostand(species25list[[i]][-c(1:2)], method = "total",1), method = "bray")
                mds.ab  <- monoMDS(stepacross(dist.ab), y = cmdscale(dist.ab, k=2),k = 2, model = "global", threshold = 0.8, maxit = 200, 
                                   weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
                temp <- cbind(species25list[[i]][c(1:2)],scores(mds.ab)[,1:2])
                listnmdsab[[i]] <- temp
                
                dist.pa <- vegdist(decostand(species25list[[i]][-c(1:2)], method = "pa",1), method = "bray")
                mds.pa  <- monoMDS(stepacross(dist.pa), y = cmdscale(dist.pa, k=2),k = 2, model = "global", threshold = 0.8, maxit = 200, 
                                   weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
                temp <- cbind(species25list[[i]][c(1:2)],scores(mds.pa)[,1:2])
                listnmdspa[[i]]<-  temp
              }
             

              
              #join ENVI and List of NMDS
              
              
              listNMDSenviab <- list()
              listNMDSenvipa <- list()
              
              for (i in seq(length(listnmdsab))){
                listNMDSenviab[[i]] <- join(listnmdsab[[i]], envi, by = c("TrNumber","subunit"), type = "left", match = "all")
                listNMDSenvipa[[i]] <- join(listnmdspa[[i]], envi, by = c("TrNumber","subunit"), type = "left", match = "all")
              }

              listNMDSenviab[[1]]$group <- rep("Pteridophyte", length( listNMDSenviab[[1]][,1]))
              listNMDSenviab[[2]]$group <- rep("Zingiberales" , length( listNMDSenviab[[2]][,1]))
              listNMDSenviab[[3]]$group <- rep("Palms"    , length( listNMDSenviab[[3]][,1]))
              listNMDSenviab[[4]]$group <- rep("Melastomataceae", length(listNMDSenviab[[4]][,1]))
              
              listNMDSenvipa[[1]]$group <- rep("Pteridophyte", length( listNMDSenvipa[[1]][,1]))
              listNMDSenvipa[[2]]$group <- rep("Zingiberales" , length( listNMDSenvipa[[2]][,1]))
              listNMDSenvipa[[3]]$group <- rep("Palms"    , length( listNMDSenvipa[[3]][,1]))
              listNMDSenvipa[[4]]$group <- rep("Melastomataceae", length(listNMDSenvipa[[4]][,1]))
              
             
              #put the dataframes together instead of having as a list
              
              unlistlistNMDSenviab <- dplyr::bind_rows(listNMDSenviab[[1]],listNMDSenviab[[2]],listNMDSenviab[[3]],listNMDSenviab[[4]])
              unlistlistNMDSenvipa <- dplyr::bind_rows(listNMDSenvipa[[1]],listNMDSenvipa[[2]],listNMDSenvipa[[3]],listNMDSenvipa[[4]])
              
              
              #prepare ENVI for join here - the name of the variables will appear in the graph of the regression trees
              #LISTS for the REGRESSION TREES
              envi2 <- envi
              colnames(envi2) <- c("TrNumber" ,"subunit","LOI","pH","Ca","K","Mg","Na","sum_of_basis","Al","P","surface","topo","lon","lat","HAND",
                                   "area_c","Slope","Drainage","max_drain","SCC","solostopoint","srtm_alt","Truni")
              
              listNMDSenviab_t <- list()
              listNMDSenvipa_t <- list()
              
              for (i in seq(length(listnmdsab))){
                listNMDSenviab_t[[i]] <- join(listnmdsab[[i]], envi2, by = c("TrNumber","subunit"), type = "left", match = "all")
                listNMDSenvipa_t[[i]] <- join(listnmdspa[[i]], envi2, by = c("TrNumber","subunit"), type = "left", match = "all")
              }
              
              listNMDSenviab_t[[1]]$group <- rep("Pteridophyte", length(listNMDSenviab_t[[1]][,1]))
              listNMDSenviab_t[[2]]$group <- rep("Zingiberales" , length(listNMDSenviab_t[[2]][,1]))
              listNMDSenviab_t[[3]]$group <- rep("Palms"    , length(listNMDSenviab_t[[3]][,1]))
              listNMDSenviab_t[[4]]$group <- rep("Melastomataceae", length(listNMDSenviab_t[[4]][,1]))
              
              listNMDSenvipa_t[[1]]$group <- rep("Pteridophyte", length(listNMDSenvipa_t[[1]][,1]))
              listNMDSenvipa_t[[2]]$group <- rep("Zingiberales" , length(listNMDSenvipa_t[[2]][,1]))
              listNMDSenvipa_t[[3]]$group <- rep("Palms"    , length(listNMDSenvipa_t[[3]][,1]))
              listNMDSenvipa_t[[4]]$group <- rep("Melastomataceae", length(listNMDSenvipa_t[[4]][,1]))
              
          

#end of the preparation of the list of NMDS and envi
              
################################################################### graphics nmds x hand #####

              #nmds vs hand
              
              nmdsHAND<- ggplot(unlistlistNMDSenvipa,aes(hand,MDS2))+
                geom_point(aes(color=surface,size=1))+
                #geom_smooth(aes(hand,MDS1,color=surface),method = "lm",se=FALSE)+
                facet_grid(group~surface,scales="free") +
                theme(strip.text.x = element_text(size=20))+
                theme(strip.text = element_text(size=25))+
                theme(legend.position="none")+
                #geom_text(aes(label=Truni))+
                theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
                ylab("NMDS1")+ 
                xlab("HAND")+
                theme(axis.title.y = element_text(size = rel(1.5)))+
                theme(axis.title.x = element_text(size = rel(1.5)))
              nmdsHAND
              
              nmdsSoil<- ggplot(unlistlistNMDSenvipa,aes(soloslinearint,MDS1))+
                geom_point(aes(color=surface,size=1))+
                #geom_smooth(aes(hand,MDS1,color=surface),method = "lm",se=FALSE)+
                facet_grid(group~surface,scales="free") +
                theme(strip.text.x = element_text(size=20))+
                theme(strip.text = element_text(size=25))+
                theme(legend.position="none")+
                #geom_text(aes(label=Truni))+
                theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
                ylab("NMDS1")+ 
                xlab("Soil Cation Concentration")+
                theme(axis.title.y = element_text(size = rel(1.5)))+
                theme(axis.title.x = element_text(size = rel(1.5)))+
                theme(axis.text = element_text(size = rel(1.5)))+
                scale_x_log10(breaks=c(0.10,0.25,1,10,50))
              nmdsSoil
              
              
              mfrow=c(1,1)
              tiff(filename="RESU_nmds_hand.tiff", height = 12*mfrow[1], width = 14*mfrow[1],units = "in",res = 300) # set the pdf file
              par(mfrow=mfrow, mar=c(0.2,0.2,0.2,0.2), oma=c(2,1,.5,0.5), mgp=c(1.7,0.6,0))
              nmdsHAND
              dev.off()
              
              mfrow=c(1,1)
              tiff(filename="RESU_nmds_soil2.tiff", height = 12*mfrow[1], width = 14*mfrow[1],units = "in",res = 300) # set the pdf file
              par(mfrow=mfrow, mar=c(0.2,0.2,0.2,0.2), oma=c(2,1,.5,0.5), mgp=c(1.7,0.6,0))
              nmdsSoil
              dev.off()

              
# end of graphics part 1

##################################################################  mixed models #####
############### models ####
              
              listNMDSmodelsPA <- list()
              listNMDSmodelsPA2 <- list()
              listAnovas <- list()
              listAnovas2 <- list()
              listhidro <- list() #this one is to calculate the number os transects and subunits used in each nmds models analysis
              
              for (i in seq(listNMDSenvipa)){
                
                hidro <- listNMDSenvipa[[i]][which(listNMDSenvipa[[i]][,"hand"]!="NA"),]
                hidro$hand <- log(1+hidro$hand)
                hidro$soloslinearint <- log(hidro$soloslinearint)
                toscale <- c("hand","decli","drain","soloslinearint", "solostopoint","srtm_alt")
                hidro[,toscale] <- scale(hidro[,toscale])
                
                listhidro[[i]] <- hidro
                
                modelNMDS <- lmer(MDS1~hand+decli+drain+soloslinearint+(1|TrNumber),data=hidro)
                listNMDSmodelsPA[[i]] <- modelNMDS 
                listAnovas[[i]] <- as.data.frame(Anova(modelNMDS)[,c(1,3)])
                listAnovas[[i]] <- round(listAnovas[[i]],2)
                
                modelNMDS2 <- lmer(MDS2~hand+decli+drain+soloslinearint+(1|TrNumber),data=hidro)
                listNMDSmodelsPA2[[i]] <- modelNMDS2 
                listAnovas2[[i]] <- as.data.frame(Anova(modelNMDS2)[,c(1,3)])
                listAnovas2[[i]] <- round(listAnovas2[[i]],2)
                
              }
              
              listaanovasNMDSmodels <- ldply(listAnovas,data.frame) 
              colnames(listaanovasNMDSmodels) <- c("Chisq", "p-value")
              
              listaanovasNMDSmodels2 <- ldply(listAnovas2,data.frame) 
              colnames(listaanovasNMDSmodels2) <- c("Chisq", "p-value")
              
              #make an output table
              groups_names <- c(rep("Ferns",4),rep("Zingiberales",4),rep("Arecaceae",4),rep("Melastomataceae",4))
              variables_names <-rep(c("HAND","Slope","Drainage","Cations"),4)
              
              Aics <- c(rep(round(AIC(listNMDSmodelsPA[[1]]),0),4),
                        rep(round(AIC(listNMDSmodelsPA[[2]]),0),4),
                        rep(round(AIC(listNMDSmodelsPA[[3]]),0),4),
                        rep(round(AIC(listNMDSmodelsPA[[4]]),0),4))
              listaanovasNMDSmodels <- cbind(variables_names,listaanovasNMDSmodels,groups_names,Aics)
              
              Aics2 <- c(rep(round(AIC(listNMDSmodelsPA2[[1]]),0),4),
                        rep(round(AIC(listNMDSmodelsPA2[[2]]),0),4),
                        rep(round(AIC(listNMDSmodelsPA2[[3]]),0),4),
                        rep(round(AIC(listNMDSmodelsPA2[[4]]),0),4))
              listaanovasNMDSmodels2 <- cbind(variables_names,listaanovasNMDSmodels2,groups_names,Aics2)
            
              write.csv(listaanovasNMDSmodels,"RESU_LMM_NMDS1.csv",row.names = TRUE)
              write.csv(listaanovasNMDSmodels2,"RESU_LMM_NMDS2.csv",row.names = TRUE)
              
              
              #how many subunits?
              lapply(listhidro,dim)
              
              #how many transects?
              length(unique(listhidro[[1]]$TrNumber))
              length(unique(listhidro[[2]]$TrNumber))
              length(unique(listhidro[[3]]$TrNumber))
              length(unique(listhidro[[4]]$TrNumber))
              
####################################################### regression tress ####

              
              
              # regression tree
              

              listTrees <- list()
              listTreescp <- list()
              
              for (j in seq(listNMDSenvipa_t)){
              for (i in unique(listNMDSenvipa_t[[j]]$surface)){
                temp <- subset(listNMDSenvipa_t[[j]], listNMDSenvipa_t[[j]]$surface==i)
                temp <- temp[which(temp$SCC!="NA"),]
                temp <- temp[which(temp$HAND!="NA"),]
                tree <- mvpart(MDS1+MDS2~SCC+HAND+Slope+Drainage,data=temp,plot.add = FALSE,cp=0.045,xv="min",margin=0.1,
                               pretty = FALSE,digits = 2,bars = FALSE,legend = FALSE)
                #nome <- unique(listNMDSenvipa_t[[j]]$group)
                #if(j==1){mtext(i,side=3,cex=1.5)}
                #if(i=="River terraces"){mtext(nome,side=2,cex=1.5)}
                listTrees[[paste(i,j)]] <- tree
                listTreescp[[paste(i,j)]] <-  tree$cptable[,1:3]
              }
              }
              
              
              mfrow=c(4,3)
              tiff(filename="Trees.tiff", height = 2.5*mfrow[1], width = 3*mfrow[1],units = "in",res = 96) # set the pdf file
              par(mfrow=mfrow, mar=c(0.1,3,4,0.3), oma=c(2,3,5,0.5), mgp=c(1.7,0.6,0))
              
              for(i in seq(listTrees)){
              rpart.plot(listTrees[[i]],type=3,cex=1.1,extra=1,box.palette = 0)
              mtext(paste("River terraces","                             ","Içá formation",
                          "                       ",
                          "Solimões formation"),side=3,outer=TRUE,line=1,cex=1.5)
              mtext(paste("Melastomataceae","           ","Zingiberales","            ","Palms","            ","Pteridophyte"),side=2,outer=TRUE,cex=1.5)
              }
              dev.off()
              

              ##teste

              listTrees <- list()
              listTreescp <- list()
              
              for (j in seq(listpatotal)){
                for (i in unique(listNMDSenvipa[[j]]$surface)){
                  temp <- subset(listNMDSenvipa[[j]], listNMDSenvipa[[j]]$surface==i)
                  temp <- temp[which(temp$soloslinearint!="NA"),]
                  temp <- temp[which(temp$hand!="NA"),]
                  tree <- mvpart(listpatotal[[j]]~log(soloslinearint)+log(1+hand)+decli,data=temp,plot.add = FALSE,cp=0.01,xv="min")
                  listTrees[[paste(i,j)]] <- tree
                  listTreescp[[paste(i,j)]] <-  tree$cptable[,1:3]
                }
              }
            
              
              plot(listNMDSenvipa[[4]]$MDS1,log(listNMDSenvipa[[4]]$soloslinearint))
              
########################################################################### graphic species tolerances along the HAND gradient ####          
              
              coefresult <- list()
              
              for(i in seq(length(listrelaabutotal))){
                temp <- which(listspENVI[[i]][,"hand"]!="NA")
                noColzero <- which(colSums(listrelaabutotal[[i]][temp,])!=0)
                
                fit_25 <- WA(listrelaabutotal[[i]][temp,noColzero], listspENVI[[i]][temp,"hand"], tolDW=TRUE) 
                coefresult[[i]] <- as.data.frame(coefficients(fit_25))
                coefresult[[i]]$tol_min <- coefresult[[i]]$Optima-coefresult[[i]]$Tolerances
                coefresult[[i]]$tol_max <- coefresult[[i]]$Optima+coefresult[[i]]$Tolerances 
                coefresult[[i]]$group   <- names(species25list)[[i]]
                coefresult[[i]]$species <- rownames(coefresult[[i]])
                
                temp2 <- coefresult[[i]]$tol_min
                temp2[temp2<0] <- 0
                coefresult[[i]]$tol_min <- temp2
              }
              
              unlistcoefresult <- plyr::ldply(coefresult) #transform the list into a single table
              
              
              coefgraph_f <- ggplot(coefresult[[1]],aes(x=reorder(species,Optima),y=Optima,ymax=tol_max,ymin=tol_min)) +
                geom_pointrange(position=position_dodge(width=0.5),size=0.8) +
                ylab("")+
                xlab("")+
                theme_bw() +
                coord_flip()+
                theme(axis.text.y = element_blank())+
                theme(axis.title = element_text(size = rel(1.8)),panel.grid.major = element_blank(),panel.grid.minor=element_blank())+
                theme(legend.position="none")+
                facet_wrap(~group,scales="free_x")+
                theme(strip.text.x = element_text(size=12, face="bold"))
              coefgraph_f
              
              coefgraph_z <- ggplot(coefresult[[2]],aes(x=reorder(species,Optima),y=Optima,ymax=tol_max,ymin=tol_min)) +
                geom_pointrange(position=position_dodge(width=0.5),size=0.8) +
                ylab("Soil Moisture (HAND)")+
                xlab("")+
                theme_bw() +
                coord_flip()+
                theme(axis.text.y = element_blank())+
                theme(axis.title = element_text(size = rel(1.8)),panel.grid.major = element_blank(),panel.grid.minor=element_blank())+
                theme(legend.position="none")+
                facet_wrap(~group,scales="free_x")+
                theme(strip.text.x = element_text(size=12, face="bold"))
              coefgraph_z
              
              coefgraph_p <- ggplot(coefresult[[3]],aes(x=reorder(species,Optima),y=Optima,ymax=tol_max,ymin=tol_min)) +
                geom_pointrange(position=position_dodge(width=0.5),size=0.8) +
                ylab("")+
                xlab("")+
                theme_bw() +
                coord_flip()+
                theme(axis.text.y = element_blank())+
                theme(axis.title = element_text(size = rel(1.8)),panel.grid.major = element_blank(),panel.grid.minor=element_blank())+
                theme(legend.position="none")+
                facet_wrap(~group,scales="free_x")+
                theme(strip.text.x = element_text(size=12, face="bold"))
              coefgraph_p
              
              mfrow=c(1,3)
              tiff(filename="RESU_figure_species_tolerances.tiff", height = 7.5*mfrow[1], width = 17*mfrow[1],units = "in",res = 600) # set the pdf file
              par(mfrow=mfrow, mar=c(0.2,0.2,0.2,0.2), oma=c(2,1,.5,0.5), mgp=c(1.7,0.6,0))
              multiplot(coefgraph_f,coefgraph_z,coefgraph_p,cols = 3)
              dev.off()
              
              #end of the tolerances analysis
              
       
              