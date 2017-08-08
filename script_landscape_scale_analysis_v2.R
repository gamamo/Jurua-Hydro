#setwd("C:/workspace gabriel/hidrologia do jurua/camilo R script june 2017")
setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/camilo R script june 2017")


#Load the bilinear extractions made in R by Camilo
      hand53_ream_m0 <-  read.csv("transect_adjust_m0_ream_hand.csv",sep = ";")  
        head(hand53_ream_m0)
        hand53_ream_m0$sub <- rep(seq(1:20),length(unique(hand53_ream_m0$tr)))
      hand53_ream_01<-  read.csv("transect_adjust_01_ream_hand.csv",sep = ";")  
        head( hand53_ream_01)


getwd()
#setwd("C:/workspace gabriel/hidrologia do jurua/Analyses")
setwd("C:/workspace_Gabriel_Moulatlet/Hidrologia do jurua/Analyses")


#usar os transectos encontrados em hand53_ream_01$tr menos o tr 784
#os tr 743,742,745,764,770,775,785 tiveram péssimos ajustes e no momento tb serão excluídos
        toexclude <- as.factor(c(743,742,745,764,770,775,784,785,805))
        trvector <- as.factor(unique(hand53_ream_01$tr))
        trvector <- trvector[!trvector %in% toexclude]
        trvector  <- as.numeric(as.vector(unique(trvector)))

#rodar o nmds for all the transects

  
  fern25_sel <- fern25[fern25$TrNumber %in% trvector,]
  
    #exclude weird subunits
    id <- seq(1:length(fern25_sel[,1]))
    ex <- cbind(fern25_sel,id)
    ex1 <- ex[ex$TrNumber=="777" & ex$subunit=="15","id"] # bizarra no mds
    ex2 <- ex[ex$TrNumber=="798" & ex$subunit=="12","id"] # bizarra no mds
    ex3 <- ex[ex$TrNumber=="766" & ex$subunit==c("1","2"),"id"] #outliers de hand
    ex4 <- ex[ex$TrNumber=="788" & ex$subunit=="1","id"]  #outliers de hand
    ex5 <- ex[ex$TrNumber=="802" & ex$subunit=="19","id"] #outliers de hand
    exx <- c(ex1,ex2,ex3,ex4,ex5)
    
  
  dist.ab <- vegdist(decostand(fern25_sel[-exx,-c(1:2)], method = "total",1), method = "bray")
  mds.ab  <- monoMDS(dist.ab, y = cmdscale(dist.ab, k=1),k = 1, model = "global", threshold = 0.8, maxit = 200, 
                     weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999)
            dist.final <- vegdist(scores(mds.ab),"euclid")
            summary(lm(dist.ab~dist.final))
            #plot(mds.ab)
            #scores(mds.ab)
            
  dist.pa <- vegdist(decostand(fern25_sel[-exx,-c(1:2)], method = "pa",1), method = "bray")
  mds.pa  <- monoMDS(dist.pa, y = cmdscale(dist.pa, k=1),k = 1, model = "global", threshold = 0.8, maxit = 200, 
                     weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999)
            dist.final <- vegdist(scores(mds.pa),"euclid")
            summary(lm(dist.pa~dist.final))
            #plot(mds.pa)
            #scores(mds.pa)
            
  
  
  #prepare the fern_mds dataset
    #get the geological info
  geoID_sel <- geoID[geoID$transect %in% trvector,] 
  geoID_sel$megageo <- NA
  geoID_sel[which(geoID_sel$surface=="Hills"),"megageo"]   <- "Ica"
  geoID_sel[which(geoID_sel$surface=="Terrace"),"megageo"] <- "Ica"
  geoID_sel[which(geoID_sel$surface=="Pebas"),"megageo"]   <- "Solimoes"
  
  #assign it to the fern_mds scores
  fern_mds <- cbind(fern25_sel[-exx,c(1:2)],scores(mds.pa))
  length(scores(mds.pa)[,1])
  
  for(i in unique(fern_mds$TrNumber)){
    fern_mds[fern_mds$TrNumber==i, "surface"] <- geoID_sel[geoID_sel$transect==i, "surface"]
    fern_mds[fern_mds$TrNumber==i, "megageo"] <- geoID_sel[geoID_sel$transect==i, "megageo"] 
  }
  
  fern_mds$Truni <- as.numeric(paste(fern_mds$TrNumber,fern_mds$subunit,sep = ""))
  
  #########################################
  
  #load the topography data
  topography <- read.csv("topography.csv", header = TRUE, sep = ",", dec = ".", na.strings = NA)
  head(topography)
  
  dist <- c(1:20 * 25)
  
  #interpolate the topography
  listtopointerpolation <- list()
  for (k in unique(trvector)){
    topo_sel <- subset(topography, topography$tr==k)
    topoapprox <- approx(topo_sel$sub,topo_sel$topo,rev(dist))
    
    matrixtopo <- matrix(NA, ncol=3, nrow = length(topoapprox$y))
    matrixtopo[,1] <- topoapprox$y
    matrixtopo[,2] <- rep(k, times=length(topoapprox$y))
    matrixtopo[,3] <- rep(1:20)
    numerator <- which(trvector==k)
    listtopointerpolation[[numerator]] <- matrixtopo
  }
  topo_sel <- ldply(listtopointerpolation,data.frame)
  colnames(topo_sel) <- c("topo","trnumber","sub")
    topo_sel$Truni <- as.numeric(paste(topo_sel$trnumber,topo_sel$sub,sep = ""))
    
    topo_sel <- topo_sel[topo_sel$Truni %in% fern_mds$Truni,]
    head(topo_sel)
    
    #calculate the differences of topo
  
  #############################################################################
  

  #load hands
  

  hand53_ream_m0_sel <- hand53_ream_m0[hand53_ream_m0$tr %in% trvector,]
  
            for(i in unique(hand53_ream_m0_sel$tr)){
              hand53_ream_m0_sel[hand53_ream_m0_sel$tr==i, "surface"] <- geoID_sel[geoID_sel$transect==i, "surface"] 
              hand53_ream_m0_sel[hand53_ream_m0_sel$tr==i, "megageo"] <- geoID_sel[geoID_sel$transect==i, "megageo"] 
            }
  
          hand53_ream_m0_sel$Truni <- as.numeric(paste(hand53_ream_m0_sel$tr,hand53_ream_m0_sel$sub,sep = ""))
    
    hand53_ream_m0_sel <- hand53_ream_m0_sel[hand53_ream_m0_sel$Truni %in% fern_mds$Truni,]

  if(F){      
  hand53_ream_01_sel <- hand53_ream_01[hand53_ream_01$tr %in% trvector, ]
            for(i in unique(hand53_ream_01_sel$tr)){
              hand53_ream_01_sel[hand53_ream_01_sel$tr==i, "surface"] <- geoID_sel[geoID_sel$transect==i, "surface"] 
              hand53_ream_01_sel[hand53_ream_01_sel$tr==i, "megageo"] <- geoID_sel[geoID_sel$transect==i, "megageo"] 
            }
    hand53_ream_01_sel <- hand53_ream_01_sel[-exx,]
  }
    
    
    #to plot
    toplot <- cbind(fern_mds,hand53_ream_m0_sel,topo_sel$topo)
      colnames(toplot)[21] <- "topo"
      
    toplot$TrNumber <- as.factor(toplot$TrNumber)
    head(toplot)
    toplot <- toplot[,which(duplicated(t(toplot))==FALSE)]
    
    # assign the topographic differences of each transect to a ranking
    toplot$topodifford <- NA 
    
    sel_list <- list()
    for(i in unique(toplot$TrNumber)){
      sel <- subset(toplot,toplot$TrNumber==i)
      sel$topodifford <- rep(dici_topodiff[dici_topodiff$TrNumber==i,"topodifford"],times=length(sel[,1]))
      #numerator <- which(toplot$TrNumber==i)
      sel_list [[i]] <- sel
    }
    toplot <- ldply(sel_list,data.frame)
      toplot <- toplot[,-1]
      toplot <- toplot[order(toplot$topodifford),]
      
      
      toplot$topodifford <- as.factor(toplot$topodifford )
      #toplot$drain       <- as.factor(toplot$drain)
      head(toplot)
      
      subset(toplot,toplot$topodifford==14)
      

    #plots
      g <- ggplot(toplot,aes(hand,MDS1))+
          geom_point(aes(size=3,color=surface))+
          geom_smooth(aes(hand,MDS1),method = "loess",se=FALSE)+
          facet_wrap(~topodifford, scales="free") +
          theme(strip.text.x = element_text(size=20))+
          theme(strip.text = element_text(size=25))
      g
      g1<- ggplot(toplot,aes(topo,MDS1))+
        geom_point(aes(color=surface,size=3))+
        geom_smooth(aes(topo,MDS1,color=surface),method = "lm",se=FALSE)
      g1
      g2<- ggplot(toplot,aes(hand,decli))+
        geom_point(aes(color=surface,size=3))+
        geom_smooth(aes(hand,decli,color=surface),method = "lm",se=FALSE)+
        facet_wrap(~TrNumber, scales="free") +
        theme(strip.text.x = element_text(size=20))+
        theme(strip.text = element_text(size=25))+
        theme(legend.position="none")
      g2

      ##
      ## mixed models based on hydrological features only
      head(toplot)
      mf1 <- lm(MDS1~hand*decli*drain, data=toplot);summary(mf1)
        plot(mf1)
        Ef1 <- rstandard(mf1)
        boxplot(Ef1~TrNumber, data=toplot)
        abline(0,0)
      mf1.gls <- gls(MDS1~hand*decli*drain, data=toplot);summary(mf1.gls)
      mf1.lme <- lme(MDS1~hand*decli*drain, random = ~1+hand|TrNumber ,data=toplot);summary(mf1.lme)
      anova(mf1.gls,mf1.lme)
      

      #FINAL MODEL
      mf2.lme <- lme(MDS1~hand*decli*drain, random = ~1+hand|TrNumber,method = "REML" ,data=toplot);summary(mf2.lme)
      plot(allEffects(mf2.lme))

      #correlations between observations from the same tr
      0.943766^2/(0.943766^2+0.3449088^2)

      #validation
      Ef2 <- resid(mf2.lme,type="normalized")
      F2  <- fitted(mf2.lme) 
      plot(F2,Ef2)
      boxplot(Ef2~TrNumber, data=toplot)
      abline(0,0)
      
      plot(mf2.lme, megageo~resid(.),abline=c(0,0))
      plot(mf2.lme, resid(., type="p")~fitted(.)|TrNumber)
      plot(mf2.lme, MDS1 ~fitted(.)|TrNumber )
      
      
    
 ########################################################################
 ########################################################################
 
  # combine solo features and hydrological features
  #solo e hand
  solo_hand <- solo[solo$TrNumber %in% trvector,]
  
  solo_hand$Truni <- as.numeric(paste(solo_hand$TrNumber,solo_hand$subunit,sep = ""))
    exsh <- solo_hand[solo_hand$TrNumber=="798" & solo_hand$subunit=="12","Truni"]
    which(solo_hand$Truni==exsh)
    solo_hand <- solo_hand[-which(solo_hand$Truni==exsh),]
  hand53_ream_m0_sel$Truni <- as.numeric(paste(hand53_ream_m0_sel$tr,hand53_ream_m0_sel$sub,sep=""))
  fern_mds$Truni <- as.numeric(paste(fern_mds$TrNumber,fern_mds$subunit,sep=""))
  topo_sel$Truni <- as.numeric(paste(topo_sel$trnumber,topo_sel$sub,sep=""))
  
  for (i in solo_hand$Truni){
    solo_hand[solo_hand$Truni==i, "hand" ] <-  hand53_ream_m0_sel[hand53_ream_m0_sel$Truni==i, "hand"]
    solo_hand[solo_hand$Truni==i, "MDS1" ] <-  fern_mds[fern_mds$Truni==i, "MDS1"]
    solo_hand[solo_hand$Truni==i, "topo" ] <-  topo_sel[topo_sel$Truni==i, "topo"]
    solo_hand[solo_hand$Truni==i, "decli" ] <-  hand53_ream_m0_sel[hand53_ream_m0_sel$Truni==i, "decli"]
    solo_hand[solo_hand$Truni==i, "drain" ] <-  hand53_ream_m0_sel[hand53_ream_m0_sel$Truni==i, "drain"]
    solo_hand[solo_hand$Truni==i, "surface" ] <-  hand53_ream_m0_sel[hand53_ream_m0_sel$Truni==i, "surface"]
    solo_hand[solo_hand$Truni==i, "megageo" ] <-  hand53_ream_m0_sel[hand53_ream_m0_sel$Truni==i, "megageo"]
  }



#mixed models Zuur protocol
  #explore the data
  #model HAND and SOIL FERTILITY
  
    head(solo_hand)
  
    m1 <- lm(MDS1 ~ hand*sum_of_basis*decli, data = solo_hand )
    summary(m1)
    plot(m1, select=c(1))
    E <- rstandard(m1)
    boxplot(E~TrNumber,data=solo_hand)
    
    m.gls <- gls(MDS1 ~ hand+log(sum_of_basis), data = solo_hand)
    summary(m.gls)
    
    ### MODELO FINAL
    m1.lme <- lme(MDS1 ~ hand+log(sum_of_basis), random = ~1|TrNumber, data = solo_hand)
    summary(m1.lme)
    plot(allEffects(m1.lme))
    
    anova(m.gls,m1.lme)
    
    E2 <- resid(m1.lme, type="normalized")
    F2 <- fitted(m1.lme)
    plot(F2,E2)
    boxplot(E2~TrNumber, data=solo_hand)
    abline(0,0)
    plot(log(solo_hand$sum_of_basis),E2)
    plot(solo_hand$hand,E2)


  
  # regression tree
  data <- solo_hand[,c("MDS1","sum_of_basis","hand")]
  data <- decostand(data,method="standardize",MARGIN = 2)
  tree <- mvpart(MDS1~sum_of_basis,data=data)
  summary(tree)
    
 
  fern_mds$Truni %in% solo_hand$Truni
  solo_hand$Truni %in%fern_mds$Truni
  unique(fern_mds$TrNumber)
  unique(solo_hand$TrNumber)
  
  g2 <- ggplot(solo_hand,aes(sum_of_basis,MDS1))+
    geom_point(aes(size=3))+
    geom_smooth(aes(sum_of_basis,MDS1),method = "lm",se=FALSE)+
    facet_wrap(~TrNumber, scales="free") +
    theme(strip.text.x = element_text(size=20))+
    theme(strip.text = element_text(size=25))
  g2
  
  g3 <- ggplot(solo_hand,aes(hand,MDS1))+
    geom_point(aes(size=3,color=surface))+
    geom_smooth(aes(hand,MDS1),method = "lm",se=FALSE)+
    facet_wrap(~TrNumber, scales="free") +
    theme(strip.text.x = element_text(size=20))+
    theme(strip.text = element_text(size=25))
  g3
  

  
 ########################################################################################
  #and interpolate the data for each transect
  colnames(solo_hand)
  
  MDSrotate(solo_hand[,c("MDS1","MDS2")],solo_hand$topo)
  
  ordisurf(solo_hand[,c("MDS1","MDS2")],solo_hand$topo,choices = c(1,2),bs="tp",col="forestgreen",main = "soil fertility",knots = 10,cex=log(2+solo_hand$sum_of_basis),pch=19)
  

  #rotate the mds to plot the  ordisurf function
  mds.ab.rot <- MDSrotate(mds.pa,topo_sel$topo)

  
  # run the trend line
  par(mfrow=c(1,2))
  ordisurf(mds.ab.rot,topo_sel$topo,choices = c(1,2),bs="tp",col="forestgreen",main = "soil fertility",knots = ,cex=log(2+solo_sel$sum_of_bases),pch=19)
  ordisurf(mds.ab.rot,topo_sel$topo,choices = c(1,2),bs="tp",col="forestgreen",main = "soil moisture",knots = 10,
           cex=hand53_ream_m0_sel$hand/15,pch=19)
  ordisurf(mds.ab.rot,topo_sel$topo,choices = c(1,2),bs="tp",col="forestgreen",main = paste("hand 01"),knots = 10,
           cex=hand53_ream_m0_sel/20,pch=19)
  
  
  
  
  #plot trend line summary
  summary(ordisurf(x = mds.ab.rot, y = topo_sel$topo, knots = 10, plot = FALSE))
  rsq <- summary(ordisurf(x = mds.ab.rot, y = topo_sel$topo, knots = 10, plot = FALSE))$r.sq
  legend("bottomright", y=NULL, legend= paste("R2","=", round(rsq,2)),cex=2,bty="n" )
  
  dev.off()
  
  #chech if variables correlate
  plot(topo_sel$topo,hand53_ream_m0_sel)
  plot(topo_sel$topo,hand53_ream_01_sel)
  plot(topo_sel$topo,log(solo_sel$sum_of_base))
  plot(log(solo_sel$sum_of_base),hand53_ream_m0_sel)
  plot(log(solo_sel$sum_of_base),hand53_ream_01_sel)
  
  #check if the first axis correlate with the envi
  plot(topo_sel$topo,scores(mds.ab)[,1])
  plot(solo_sel$sum_of_base,scores(mds.ab)[,1])
  plot(hand53_ream_m0_sel,scores(mds.ab)[,2])
  plot(hand53_ream_01_sel,scores(mds.ab)[,1])
  plot(scores(mds.ab)[,1],scores(mds.ab)[,2],cex=hand53_ream_m0_sel/15)
  
  m <- nls(scores(mds.ab)[,1]~a*hand53_ream_m0_sel/(b+hand53_ream_m0_sel))
  cor(scores(mds.ab)[,1],predict(m))
  
  
  
 