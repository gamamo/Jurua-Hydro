      # Jurua analyses
      # this script will be used to run the analyses of the Jurua top-hydro paper
      # owner: Gabriel Massaine Moulatlet
      # contact: gamamo@utu.fi
      
      
      # load the relevant packages
      
      library(vegan)
      library(ggplot2)
      library(lme4)
      library(GGally)
      library(plyr)
      library(dplyr)
      library(nlme)
      library(effects)
      library(car)
      library(pROC)
      library(mvpart)
      library(rioja)
      library(palaeoSig)
      library(rpart.plot)
      library(rpart)
      

      
      # this script starts with the preparation of table that gets the relative elevational of the transects in each
      # geological formation, then come the graphics that relate environment and species
      
###################################### load the environmental data ####
      
      getwd()
      #setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses/dados ambientais originais")
      setwd("F:/workspace gabriel/Hidrologia do jurua/Analyses/dados ambientais originais") 
      dir()
      
      rm(list = ls())
      
      #load GEO list
      geoID <- read.csv("geoID.csv", stringsAsFactors = FALSE)
        geoID[which(geoID$surface=="Terrace"),"surface"] <- "Rivers Terraces"
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
                  topoapprox <- approx(topo_sel$sub,topo_sel$topo,rev(dist))
                  
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
      setwd("F:/workspace gabriel/Hidrologia do jurua/camilo R script june 2017")
      #setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/camilo R script june 2017")
        
        
        #Load the bilinear extractions made in R by Camilo
      hand53_ream_m0 <-  read.csv("transect_adjust_m0_ream_hand.csv",sep = ";")  
          head(hand53_ream_m0)
          hand53_ream_m0$sub <- rep(seq(1:20),length(unique(hand53_ream_m0$tr)))
      hand53_ream_01<-  read.csv("transect_adjust_01_ream_hand.csv",sep = ";")  
          head( hand53_ream_01)
        
      
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
            envi[envi$tr==i, c("lon","lat","hand","area_c","decli", "drain", "max_drain")] <- hand53_ream_m0[hand53_ream_m0$tr==i,
                                                                                                              c("lon","lat","hand",
                                                                                                                "area_c","decli", "drain", "max_drain")]
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
            terraces <- envi[which(envi$surface == "Rivers Terraces"     ),"TrNumber"]
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
                        aggregate(envi$decli, list(envi$TrNumber), mean       )[,-1 ]
                )
              
                colnames(statsminmax) <- c("transect", "min_topo" , "max_topo", "mean_topo", "min_hand", "max_hand" , "mean_hand",
                                           "min_cation", "max_cation" , "mean_cation","min_slope", "max_slope" , "mean_slope" )
                
                
                table1 <- matrix(data = NA, nrow = 3, ncol = 4) # creat the table 
                  rownames(table1) <- c("Içá formation" , "Solimões formation", "Rivers Terraces")
                  colnames(table1) <- c("Topography" , "Soil Moisture (HAND)", "Soil Fertility (Cation concentration)", "Slope"    )
                 
                  
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
                  }
                  print(table1)
                  
                  #setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses")
                  setwd("F:/workspace gabriel/Hidrologia do jurua/Analyses") 
                  
                  write.csv(table1, "RESUtable1.csv",row.names = TRUE)
                  
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
                  
                  #exclude weird subunits
                  id <- seq(1:length(fern25[,1]))
                  ex <- cbind(fern25,id)
                  ex1 <- ex[ex$TrNumber=="777" & ex$subunit=="15","id"] # bizarra no mds
                  ex2 <- ex[ex$TrNumber=="798" & ex$subunit=="12","id"] # bizarra no mds
                  ex3 <- ex[ex$TrNumber=="766" & ex$subunit==c("1","2"),"id"] #outliers de hand
                  ex4 <- ex[ex$TrNumber=="788" & ex$subunit=="1","id"]  #outliers de hand
                  ex5 <- ex[ex$TrNumber=="802" & ex$subunit=="19","id"] #outliers de hand
                  ex6 <- ex[ex$TrNumber=="777" & ex$subunit=="16","id"] #outliers de hand
                  exx <- c(ex1,ex2,ex3,ex4,ex5,ex6)
                  
                  fern25 <- fern25[-exx,] 
                  fern25 <- fern25[fern25$TrNumber %in% trvector,]
                
            zing25 <- read.csv("zingdata_wide_gmm_v1.csv"    , stringsAsFactors = FALSE)
                  zing25 <- zing25[zing25$TrNumber %in% trvector,]
            

            palm25 <- read.csv("palms_jurua_subunit_gmm1.csv", stringsAsFactors = FALSE) # no need to delete rows for the palms: checked before
                  palm25 <- palm25[palm25$TrNumber %in% trvector,]
                
            # creat a list if the species data             
              species25list <- list(fern25, zing25, palm25) 
              names(species25list) <- c("Ferns", "Zingiberales", "Arecaceae")
              
         
            # call the environmental data
            # each plant object has a different number of rows. It has to be adequated in the environmental data: use the join::plyr 
            # calculate species optima and toleraces to the HAND gradient
              
              listrelaabutotal <- list() #only abu
              listspENVI       <- list() #join pa and envi
              
              
              for(i in seq(length(species25list))){
                sprela <- decostand(species25list[[i]][-c(1:2)], method = "total",1)   #to calculate relative abundances
                
                listrelaabutotal[[i]] <- sprela
                listspENVI[[i]]       <- join(species25list[[i]], envi, by = c("TrNumber","subunit"), type = "left", match = "all") # and merge with envi
              }

#end of species importing and envi matching

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

###################################################################### Run ordinations e prepare the table for the analysis #####
      
      
             # run the NMDS ordinations
             # first calculate the relative abundances
             # then run the NMDS and store the first 2 axis of each plant group in a separate list
       
              listnmds <- list()
              
              for (i in seq(length(species25list))){
                
                dist.ab <- vegdist(decostand(species25list[[i]][-c(1:2)], method = "total",1), method = "bray")
                mds.ab  <- monoMDS(dist.ab, y = cmdscale(dist.ab, k=1),k = 1, model = "global", threshold = 0.8, maxit = 200, 
                                   weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
                listnmds[[i]] <- species25list[[i]][c(1:2)]
                listnmds[[i]]$MDSabu <- scores(mds.ab)[,1]
                
                dist.pa <- vegdist(decostand(species25list[[i]][-c(1:2)], method = "pa",1), method = "bray")
                mds.pa  <- monoMDS(dist.pa, y = cmdscale(dist.pa, k=1),k = 1, model = "global", threshold = 0.8, maxit = 200, 
                                   weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999) 
                listnmds[[i]]$MDSpa <- scores(mds.pa)[,1]
              }
              
              
          
              #join ENVI and List of NMDS
              listNMDSenvi <- list()
              
              for (i in seq(length(listnmds))){
                listNMDSenvi[[i]] <- join(listnmds[[i]], envi, by = c("TrNumber","subunit"), type = "left", match = "all")
              }

              listNMDSenvi[[1]]$group <- rep("Ferns", length( listNMDSenvi[[1]][,1]))
              listNMDSenvi[[2]]$group <- rep("Zingiberales" , length( listNMDSenvi[[2]][,1]))
              listNMDSenvi[[3]]$group <- rep("Arecaceae"    , length( listNMDSenvi[[3]][,1]))
              
              #put the dataframes together instead of having as a list
              
              unlistlistNMDSenvi <- dplyr::bind_rows(listNMDSenvi[[1]],listNMDSenvi[[2]],listNMDSenvi[[3]])

#end of the preparation of the list of NMDS and envi
              
################################################################### graphics nmds x hand #####
              colnames(unlistlistNMDSenvi)
              
              #nmds vs hand
              
              nmdsHAND<- ggplot(unlistlistNMDSenvi,aes(hand,MDSpa))+
                geom_point(aes(color=surface,size=1))+
                geom_smooth(aes(hand,MDSpa,color=surface),method = "lm",se=FALSE)+
                facet_wrap(group~surface,scales="free") +
                theme(strip.text.x = element_text(size=20))+
                theme(strip.text = element_text(size=25))+
                theme(legend.position="none")+
                #geom_text(aes(label=Truni))+
                theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
                ylab("NMDS")+ 
                xlab("Soil Moisture")+
                theme(axis.title.y = element_text(size = rel(1.8)))+
                theme(axis.title.x = element_text(size = rel(1.8)))
              nmdsHAND
              
              mfrow=c(1,1)
              tiff(filename="RESU_nmds_hand.tiff", height = 10*mfrow[1], width = 10*mfrow[1],units = "in",res = 300) # set the pdf file
              par(mfrow=mfrow, mar=c(0.2,0.2,0.2,0.2), oma=c(2,1,.5,0.5), mgp=c(1.7,0.6,0))
              nmdsHAND
              dev.off()

              
# end of graphics part 1

################################################################## hydrological mixed models #####
###############hidro model for ferns ####
              
              hidro <- listNMDSenvi[[1]][which(listNMDSenvi[[1]][,"hand"]!="NA"),]
              
              mf1 <- lm(MDSpa~hand*decli*drain, data=hidro);summary(mf1)
              Ef1 <- rstandard(mf1)
              boxplot(Ef1~TrNumber, data=hidro)
              abline(0,0)
              mf1.gls <- gls(MDSpa~hand*decli*drain, data=hidro);summary(mf1.gls)
              mf1.lme <- lme(MDSpa~hand*decli*drain, random = ~1+surface|TrNumber ,data=hidro);summary(mf1.lme)
              anova(mf1.gls,mf1.lme)
              
              #FINAL MODEL
              mf2.lme <- lme(MDSpa~hand*decli, random = ~1+ surface|TrNumber,method = "REML" ,data=hidro);summary(mf2.lme)
              plot(allEffects(mf2.lme),rug=FALSE)
              
              hidro_output_ferns <- as.data.frame(Anova(mf2.lme)[,c(1,3)])

              #validation
              Ef2 <- resid(mf2.lme,type="normalized")
              F2  <- fitted(mf2.lme) 
              plot(F2,Ef2)
              boxplot(Ef2~TrNumber, data=hidro)
              boxplot(Ef2~surface, data=hidro)
              abline(0,0)
              
              plot(mf2.lme, surface~resid(.),abline=c(0,0))
              plot(mf2.lme, resid(., type="p")~fitted(.)|TrNumber)
              plot(mf2.lme, MDSpa ~fitted(.)|TrNumber )
              
#################hidro model for zing ####
              
              hidro <- listNMDSenvi[[2]][which(listNMDSenvi[[2]][,"hand"]!="NA"),]
              
              mf1 <- lm(MDSpa~hand*decli*drain, data=hidro);summary(mf1)
              Ef1 <- rstandard(mf1)
              boxplot(Ef1~TrNumber, data=hidro)
              abline(0,0)
              mf1.gls <- gls(MDSpa~hand*decli*drain, data=hidro);summary(mf1.gls)
              mf1.lme <- lme(MDSpa~hand*decli*drain, random = ~1+surface|TrNumber ,data=hidro);summary(mf1.lme)
              anova(mf1.gls,mf1.lme)
              
              mf2.lme <- lme(MDSpa~hand*decli, random = ~1+ surface|TrNumber,method = "REML" ,data=hidro);summary(mf2.lme)
              plot(allEffects(mf2.lme),rug=FALSE)
              
              hidro_output_zings <- as.data.frame(Anova(mf2.lme)[,c(1,3)])
              
              #validation
              Ef2 <- resid(mf2.lme,type="normalized")
              F2  <- fitted(mf2.lme) 
              plot(F2,Ef2)
              boxplot(Ef2~TrNumber, data=hidro)
              boxplot(Ef2~surface, data=hidro)
              abline(0,0)
              
              plot(mf2.lme, surface~resid(.),abline=c(0,0))
              plot(mf2.lme, resid(., type="p")~fitted(.)|TrNumber)
              plot(mf2.lme, MDSpa ~fitted(.)|TrNumber )
              
#################hidro model for palms ####
              
              hidro <- listNMDSenvi[[3]][which(listNMDSenvi[[3]][,"hand"]!="NA"),]
              
              mf1 <- lm(MDSpa~hand*decli*drain, data=hidro);summary(mf1)
              Ef1 <- rstandard(mf1)
              boxplot(Ef1~TrNumber, data=hidro)
              abline(0,0)
              mf1.gls <- gls(MDSpa~hand*decli*drain, data=hidro);summary(mf1.gls)
              mf1.lme <- lme(MDSpa~hand*decli*drain, random = ~1+surface|TrNumber ,data=hidro);summary(mf1.lme)
              anova(mf1.gls,mf1.lme)
              
              mf2.lme <- lme(MDSpa~hand*decli, random = ~1+ surface|TrNumber,method = "REML" ,data=hidro);summary(mf2.lme)
              plot(allEffects(mf2.lme),rug=FALSE)
              
              hidro_output_palms <- as.data.frame(Anova(mf2.lme)[,c(1,3)])
              
              #validation
              Ef2 <- resid(mf2.lme,type="normalized")
              F2  <- fitted(mf2.lme) 
              plot(F2,Ef2)
              boxplot(Ef2~TrNumber, data=hidro)
              boxplot(Ef2~surface, data=hidro)
              abline(0,0)
              
              plot(mf2.lme, surface~resid(.),abline=c(0,0))
              plot(mf2.lme, resid(., type="p")~fitted(.)|TrNumber)
              plot(mf2.lme, MDSpa ~fitted(.)|TrNumber )
              
# prepare a output table
              
              outputMM1 <- cbind(hidro_output_ferns,hidro_output_zings,hidro_output_palms)

# end of the hydro modelling
              
############################################################# graphs nmds soil and hydro
              
              unlistlistNMDSenvi_sub <- unlistlistNMDSenvi[which(unlistlistNMDSenvi$sum_of_basis!="NA"),]
              unlistlistNMDSenvi_sub <- unlistlistNMDSenvi_sub[which(unlistlistNMDSenvi_sub$hand!="NA"),]

              nmdsedap1<- ggplot(unlistlistNMDSenvi_sub,aes(hand,MDSabu))+
                geom_point(aes(color=surface,size=1))+
                geom_smooth(aes(hand,MDSabu,color=surface),method = "lm",se=FALSE)+
                facet_wrap(group~surface) +
                theme(strip.text.x = element_text(size=20))+
                ylab("NMDS")+ 
                xlab("Soil Moisture")+
                theme(strip.text = element_text(size=25))+
                theme(legend.position="none")+
                #geom_text(aes(label=Truni))+
                theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
                theme(axis.title.y = element_text(size = rel(1.8)))+
                theme(axis.title.x = element_text(size = rel(1.8)))
              nmdsedap1
              
              nmdsedap2<- ggplot(unlistlistNMDSenvi_sub,aes(sum_of_basis,MDSpa))+
                geom_point(aes(color=surface,size=1))+
                geom_smooth(aes(sum_of_basis,MDSpa,color=surface),method = "lm",se=FALSE)+
                facet_wrap(group~surface) +
                theme(strip.text.x = element_text(size=20))+
                theme(strip.text = element_text(size=25))+
                theme(legend.position="none")+
                #geom_text(aes(label=Truni))+
                theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
                ylab("NMDS")+ 
                xlab("Soil Fertility")+
                theme(axis.title.y = element_text(size = rel(1.8)))+
                theme(axis.title.x = element_text(size = rel(1.8)))+
                scale_x_log10()
              nmdsedap2
              
              mfrow=c(1,1)
              tiff(filename="RESU_figure_NMDS_hand_solo.tiff", height = 17*mfrow[1], width = 17*mfrow[1],units = "in",res = 300) # set the pdf file
              par(mfrow=mfrow, mar=c(0.2,0.2,0.2,0.2), oma=c(2,1,.5,0.5), mgp=c(1.7,0.6,0))
              multiplot(nmdsedap1,nmdsedap2)
              dev.off()

#end of graphics
              ##########parei aqui!!!! fazer as regressoes apenas para as subunidades que tem solo
              
              #build a R2 table for the figure ####
              listRhand <- list()
              listRsoil <- list()
              
              for(i in seq(listNMDSenvi)){
                for(j in unique(listNMDSenvi[[1]]$surface)){
                  temp <- subset(listNMDSenvi[[i]],listNMDSenvi[[i]]$surface==j)
                  reghand <- lm(temp$MDSpa~temp$hand)
                  regsoil <- lm(temp$MDSpa~log(temp$sum_of_basis))
                  listRhand[[paste(i,j)]] <- round(summary(reghand)$adj.r.squared,3)
                  listRsoil[[paste(i,j)]] <- round(summary(regsoil)$adj.r.squared,3)
                }
              }

              unlistRhand <- plyr::ldply(listRhand)
              unlistRsoil <- plyr::ldply(listRsoil)
              
              R2edaphic <- cbind(unlistRhand,unlistRsoil$V1)
              colnames(R2edaphic)[1:3] <- c("Geological surface","adjR2 Soil Moisture","adjR2 Soil fertility")
              plantgroup <- c(rep("Ferns",3),rep("Zingiberales",3),rep("Arecaceae",3))
              R2edaphic <- cbind(plantgroup,R2edaphic)
              
              write.csv(R2edaphic,"RESU_R2_edaphic_figure.csv",row.names = FALSE)
              
############################################################## edaphic modelling mixed models
########## edaphic model ferns ####         
              
              edaph <- listNMDSenvi[[1]][which(listNMDSenvi[[1]][,"sum_of_basis"]!="NA"),]
              edaph <- edaph[which(edaph$hand!="NA"),]
              
              m1 <- lm(MDSpa ~ hand+sum_of_basis, data = edaph)
              summary(m1)
              E <- rstandard(m1)
              length(E)
              boxplot(E~edaph$TrNumber)
              
              m.gls <- gls(MDSpa ~ hand+log(sum_of_basis), data = edaph)
              summary(m.gls)
              
              ### MODELO FINAL
              m1.lme <- lme(MDSpa ~ hand*log(sum_of_basis), random = ~1|TrNumber, method = "REML",data = edaph)
              summary(m1.lme)

 
              plot(allEffects(m1.lme),rug=FALSE)
              edaph_output_ferns <- as.data.frame (Anova(m1.lme2)[,c(1,3)])

              #validation
              E2 <- resid(m1.lme, type="normalized")
              F2 <- fitted(m1.lme)
              plot(F2,E2)
              boxplot(E2~TrNumber, data=edaph)
              abline(0,0)
              
              plot(m1.lme, surface~resid(.),abline=c(0,0))
              plot(m1.lme, resid(., type="p")~fitted(.)|TrNumber)
              plot(m1.lme, MDSpa ~fitted(.)|TrNumber )
              
########## edaphic model zing #### 
              
              edaph <- listNMDSenvi[[2]][which(listNMDSenvi[[2]][,"sum_of_basis"]!="NA"),]
              edaph <- edaph[which(edaph$hand!="NA"),]
              
              m1 <- lm(MDSpa ~ hand*sum_of_basis, data = edaph)
              summary(m1)
              E <- rstandard(m1)
              length(E)
              boxplot(E~edaph$TrNumber)
              
              m.gls <- gls(MDSpa ~ hand+log(sum_of_basis), data = edaph)
              summary(m.gls)
              
              ### MODELO FINAL
              m1.lme <- lme(MDSpa ~ hand*log(sum_of_basis), random = ~1+|TrNumber, method = "REML",data = edaph)
              summary(m1.lme)
              
              plot(allEffects(m1.lme),rug=FALSE)
              edaph_output_zings <- as.data.frame (Anova(m1.lme)[,c(1,3)])
              
              #validation
              E2 <- resid(m1.lme, type="normalized")
              F2 <- fitted(m1.lme)
              plot(F2,E2)
              boxplot(E2~TrNumber, data=edaph)
              abline(0,0)
              
              plot(m1.lme, surface~resid(.),abline=c(0,0))
              plot(m1.lme, resid(., type="p")~fitted(.)|TrNumber)
              plot(m1.lme, MDSpa ~fitted(.)|TrNumber )
              
########## edaphic model palms #### 
              
              edaph <- listNMDSenvi[[3]][which(listNMDSenvi[[3]][,"sum_of_basis"]!="NA"),]
              edaph <- edaph[which(edaph$hand!="NA"),]
              
              m1 <- lm(MDSpa ~ hand*sum_of_basis, data = edaph)
              summary(m1)
              E <- rstandard(m1)
              length(E)
              boxplot(E~edaph$TrNumber)
              
              m.gls <- gls(MDSpa ~ hand+log(sum_of_basis), data = edaph)
              summary(m.gls)
              
              ### MODELO FINAL
              m1.lme <- lme(MDSpa ~ hand*log(sum_of_basis), random = ~1|TrNumber, method = "REML",data = edaph)
              summary(m1.lme)
              
              plot(log(edaph$sum_of_basis),edaph$MDSpa)
              plot(allEffects(m1.lme),rug=FALSE)
              edaph_output_palms <- as.data.frame (Anova(m1.lme)[,c(1,3)])
              
              #validation
              E2 <- resid(m1.lme, type="normalized")
              F2 <- fitted(m1.lme)
              plot(F2,E2)
              boxplot(E2~surface, data=edaph)
              abline(0,0)
              
              plot(m1.lme, surface~resid(.),abline=c(0,0))
              plot(m1.lme, resid(., type="p")~fitted(.)|TrNumber)
              plot(m1.lme, MDSpa ~fitted(.)|TrNumber )

# prepare the output table
              
              outputMM2 <- cbind(edaph_output_ferns,edaph_output_zings,edaph_output_palms)

# end of the edaphic models

####################################################### combine output tables ####
              
              outputMMfinal <- round(rbind(outputMM1,outputMM2),2)
                outputMMfinal$model <- c(rep("Hydrological",3),rep("Edaphic",3))
                groups_names <- c(rep("Ferns",2),rep("Zingiberales",2),rep("Arecaceae",2))
                outputMMfinal <- rbind(groups_names,outputMMfinal)
                
                write.csv(outputMMfinal,"RESU_outputMM.csv",row.names = TRUE)
              
####################################################### regression tress ####

              
              # regression tree

              listTrees <- list()
              listTreescp <- list()
              
              for (j in seq(listNMDSenvi)){
              for (i in unique(listNMDSenvi[[j]]$surface)){
                temp <- subset(listNMDSenvi[[j]], listNMDSenvi[[j]]$surface==i)
                temp <- temp[which(temp$sum_of_basis!="NA"),]
                temp <- temp[which(temp$hand!="NA"),]
                tree <- mvpart(MDSpa~log(sum_of_basis)+hand+decli,data=temp,plot.add = FALSE,cp=0.01,xv="min")
                listTrees[[paste(i,j)]] <- tree
                listTreescp[[paste(i,j)]] <-  tree$cptable[,1:3]
              }
              }
              names(listTrees)
           plot(temp$decli,log(temp$sum_of_basis))

           
           a<- ggplot(unlistlistNMDSenvi_sub,aes(log(sum_of_basis),decli))+
             geom_point(aes(color=surface,size=1))+
             #geom_smooth(aes(log(sum_of_basis),MDSpa,color=surface),method = "lm",se=FALSE)+
             facet_wrap(~surface,scales = "free") +
             theme(strip.text.x = element_text(size=20))+
             theme(strip.text = element_text(size=25))+
             theme(legend.position="none")+
             #geom_text(aes(label=Truni))+
             #theme(panel.background = element_rect(fill = "white", colour = "grey50"))+
             #ylab("NMDS")+ 
             #xlab("Soil Fertility")+
             #theme(axis.title.y = element_text(size = rel(1.8)))+
             theme(axis.title.x = element_text(size = rel(1.8)))
           a
              