setwd("C:/workspace gabriel/hidrologia do jurua/camilo R script june 2017")

#Load the bilinear extractions made by arcgis
if(F){
hand53_ream_m0 <- read.csv("transect_adjust_m0_reamhandARC.csv",sep=";")
head(hand53_ream_m0)
hand53_ream_01 <- read.csv("transect_adjust_01_reamhandARC.csv",sep=";")
head(hand53_ream_01)
}


#Load the bilinear extractions made in R by Camilo
hand53_ream_m0 <-  read.csv("transect_adjust_m0_ream_hand.csv",sep = ";")  
  head(hand53_ream_m0)
  hand53_ream_m0$sub <- rep(seq(1:20),length(unique(hand53_ream_m0$tr)))
hand53_ream_01<-  read.csv("transect_adjust_01_ream_hand.csv",sep = ";")  
  head( hand53_ream_01)


getwd()
setwd("C:/workspace gabriel/hidrologia do jurua/Analyses")

unique(solo$TrNumber)
length(unique(solo$TrNumber))  
unique(fern25$TrNumber) 
length(unique(fern25$TrNumber))
unique(hand53_ream_m0$tr)
length(unique(hand53_ream_m0$tr))
unique(hand53_ream_01$tr)
length(unique(hand53_ream_01$tr))
unique(topography$tr) 
length(unique(topography$tr))

#comeco do looping

#usar os transectos encontrados em hand53_ream_01$tr menos o tr 784
#os tr 743,742,745,764,770,775,785 tiveram péssimos ajustes e no momento tb serão excluídos
toexclude <- as.factor(c(743,742,745,764,770,775,784,785,805))
trvector <- as.factor(unique(hand53_ream_01$tr))
trvector <- trvector[!trvector %in% toexclude]
trvector  <- as.numeric(as.vector(unique(trvector)))

teste <- as.numeric(c(751,760,759,757,763,786,796,799,795,787,772,744,738,802))
teste <- as.numeric(c(751,760,759,757,763,786,796,799,787,772,802,755))
trvector <- trvector[trvector %in% teste]

#rodar o nmds for all the transects

#listresultsregressions <- list()
#for (k in unique(trvector)){
  
  #k=791; dev.off()
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
  mds.ab  <- monoMDS(dist.ab, y = cmdscale(dist.ab, k=2),k = 2, model = "global", threshold = 0.8, maxit = 200, 
                     weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999)
            dist.final <- vegdist(scores(mds.ab),"euclid")
            summary(lm(dist.ab~dist.final))
            plot(mds.ab)
            
  dist.pa <- vegdist(decostand(fern25_sel[-exx,-c(1:2)], method = "pa",1), method = "bray")
  mds.pa  <- monoMDS(dist.pa, y = cmdscale(dist.pa, k=2),k = 2, model = "global", threshold = 0.8, maxit = 200, 
                     weakties = TRUE, stress = 1, scaling = TRUE, pc = TRUE, smin = 1e-4, sfgrmin = 1e-7, sratmax=0.99999)
            dist.final <- vegdist(scores(mds.pa),"euclid")
            summary(lm(dist.pa~dist.final))
            plot(mds.pa)
            
  
  
  #prepare the mds dataset
  geoID_sel <- geoID[geoID$transect %in% trvector,] 
  geoID_sel$megageo <- NA
  geoID_sel[which(geoID_sel$surface=="Hills"),"megageo"]   <- "Ica"
  geoID_sel[which(geoID_sel$surface=="Terrace"),"megageo"] <- "Ica"
  geoID_sel[which(geoID_sel$surface=="Pebas"),"megageo"]   <- "Solimoes"
  
  
  fern_mds <- cbind(fern25_sel[-exx,c(1:2)],scores(mds.pa))
  plot(fern_mds$MDS2,fern_mds$MDS1)
  
  for(i in unique(fern_mds$TrNumber)){
    fern_mds[fern_mds$TrNumber==i, "surface"] <- geoID_sel[geoID_sel$transect==i, "surface"]
    fern_mds[fern_mds$TrNumber==i, "megageo"] <- geoID_sel[geoID_sel$transect==i, "megageo"] 
  }
  
    fern_mds_sel <- subset(fern_mds,fern_mds$TrNumber==k)
  

  dist <- c(1:20 * 25)
  unique(fern_mds$TrNumber)
  ##########################################
  #interpolacao de solo
  
  solo
  listsolointerpolation <- list()
  for (k in unique(trvector)){
  solo_sel <- subset(solo, solo$TrNumber==k)
  soloapprox <- approx(solo_sel$subunit*25,solo_sel$sum_of_basis,dist,method = "linear",rule=2)
  #soloapprox <- spline(solo_sel$subunit*25,solo_sel$sum_of_basis,xout=dist,method = "natural")
    matrixsolo <- matrix(NA, ncol=2, nrow = length(soloapprox$y))
    matrixsolo[,1] <- soloapprox$y
    matrixsolo[,2] <- rep(k, times=length(soloapprox$y))
    numerator <- which(trvector==k)
  listsolointerpolation[[numerator]] <- matrixsolo
  }
  unlistsolos <- ldply(listsolointerpolation,data.frame)
  
  solo_sel <- unlistsolos[-exx,]
  colnames(solo_sel) <- c("sum_of_bases","trnumber")
  summary(solo_sel)
  length(solo_sel[,1])
  
  
  #########################################
  #load hands
  

  hand53_ream_m0_sel <- hand53_ream_m0[hand53_ream_m0$tr %in% trvector,]
            for(i in unique(hand53_ream_m0_sel$tr)){
              hand53_ream_m0_sel[hand53_ream_m0_sel$tr==i, "surface"] <- geoID_sel[geoID_sel$transect==i, "surface"] 
              hand53_ream_m0_sel[hand53_ream_m0_sel$tr==i, "megageo"] <- geoID_sel[geoID_sel$transect==i, "megageo"] 
            }
    hand53_ream_m0_sel <- hand53_ream_m0_sel[-exx,]
      
  hand53_ream_01_sel <- hand53_ream_01[hand53_ream_01$tr %in% trvector, ]
            for(i in unique(hand53_ream_01_sel$tr)){
              hand53_ream_01_sel[hand53_ream_01_sel$tr==i, "surface"] <- geoID_sel[geoID_sel$transect==i, "surface"] 
              hand53_ream_01_sel[hand53_ream_01_sel$tr==i, "megageo"] <- geoID_sel[geoID_sel$transect==i, "megageo"] 
            }
    hand53_ream_01_sel <- hand53_ream_01_sel[-exx,]
    

    
    #general plot
    col = as.numeric(as.factor(fern_mds$surface))
    col = as.numeric(as.factor(fern_mds$megageo))
    plot(hand53_ream_m0_sel$hand,fern_mds$MDS1,col=col,pch=19)
    text(hand53_ream_m0_sel$hand,fern_mds$MDS1,fern_mds$subunit)
    plot(hand53_ream_01_sel$hand,fern_mds$MDS1,col=col)
    plot(fern_mds$MDS2,fern_mds$MDS1,pch=19,xlim=c(-2,2),cex=log(2+solo_sel$sum_of_bases))
    plot(fern_mds$MDS2,fern_mds$MDS1,pch=19,xlim=c(-2,2),cex=hand53_ream_01_sel$hand/10)
    plot(log(2+solo_sel$sum_of_bases),hand53_ream_01_sel$hand/10)
    plot(hand53_ream_01_sel$area_c/10^4,fern_mds$MDS1)
    plot(hand53_ream_01_sel$max_drain,fern_mds$MDS1)
    plot(hand53_ream_01_sel$max_drain,hand53_ream_01_sel$drain)
    plot(hand53_ream_01_sel$drain,hand53_ream_01_sel$decli)
    
    plot(hand53_ream_m0_sel$hand,fern_mds$MDS1,col=col,pch=19,cex=log(hand53_ream_m0_sel$area_c/10^5))
    plot(hand53_ream_m0_sel$hand,fern_mds$MDS1,col=col,pch=19,cex=hand53_ream_m0_sel$drain)
    
    reghidro <- lm(fern_mds$MDS1+fern_mds$MDS1~hand53_ream_m0_sel$hand+hand53_ream_m0_sel$drain+hand53_ream_m0_sel$decli+
                     hand53_ream_m0_sel$area_c);summary(reghidro)
    
    #to plot
    toplot <- cbind(fern_mds,hand53_ream_m0_sel)
    toplot$TrNumber <- as.factor(toplot$TrNumber)
      head(toplot)
      toplot$TrNumber
      
      g <- ggplot(toplot,aes(hand,MDS1))+
          geom_point(aes(size=3,color=TrNumber))+
          geom_smooth(aes(hand,MDS1,color=TrNumber),method = "lm",se=FALSE)+
          facet_wrap(~TrNumber, scales="free") +
          theme(strip.text.x = element_text(size=20))+
          theme(strip.text = element_text(size=25))
      g
      g1<- ggplot(toplot,aes(hand,MDS1))+
        geom_point(aes(color=TrNumber,size=3))+
        geom_smooth(aes(hand,MDS1,color=TrNumber),method = "lm",se=FALSE)
      g1
    
    #define subset criteria
      
      
      seg=10
      suf="Hills"
      suf="Pebas"
      suf="Terrace"
      suf=c("Hills","Terrace")
      geo="Ica"
      geo="Solimoes"
      k=804
      
      crit=which(hand53_ream_m0_sel$tr==k)
      
      crit=which(hand53_ream_m0_sel$hand>seg & hand53_ream_m0_sel$drain>2)
      crit=which(hand53_ream_m0_sel$drain>0 )
      crit=which(hand53_ream_m0_sel$hand<seg)
      crit=which(hand53_ream_m0_sel$surface==suf &hand53_ream_m0_sel$hand>seg & hand53_ream_m0_sel$drain>1)
      crit=which(hand53_ream_m0_sel$megageo==geo)

      #prepare the objects
      solo_sel2 <- solo_sel[crit,]
        
      hand53_ream_m0_sel_sub <- subset(hand53_ream_m0_sel, hand53_ream_m0_sel$tr==k)
      hand53_ream_01_sel_sub <- subset(hand53_ream_01_sel, hand53_ream_01_sel$tr==k)
      
        hand53_ream_m0_sel_seg <- hand53_ream_m0_sel[crit,]
        hand53_ream_01_sel_seg <- hand53_ream_01_sel[crit,]
        

      #plot drain
    
      plot(hand53_ream_m0_sel_seg$hand,fern_mds$MDS1[crit])  
      col=hand53_ream_m0_sel_seg[which(hand53_ream_m0_sel_seg$drain==2),"hand"]
      plot(hand53_ream_m0_sel_seg$hand,fern_mds$MDS1[crit],cex=2,col=hand53_ream_m0_sel_seg$drain,pch=19)
      
      text(hand53_ream_m0_sel_seg$hand,fern_mds$MDS1[crit],hand53_ream_m0_sel_seg$drain)
      text(hand53_ream_m0_sel_seg$hand,fern_mds$MDS1[crit],fern_mds[crit,]$TrNumber)
      abline(lm(fern_mds$MDS1[crit]~hand53_ream_m0_sel_seg$hand))
      plot(hand53_ream_01_sel_seg$hand,fern_mds$MDS1[crit],cex=2,col=hand53_ream_01_sel_seg$drain,pch=19)

      
      #plot max order
      plot(hand53_ream_m0_sel_seg$hand,fern_mds$MDS1[crit],col=hand53_ream_m0_sel_seg$max_drain,pch=19)
      plot(hand53_ream_01_sel_seg$hand,fern_mds$MDS1[crit],cex=hand53_ream_01_sel_seg$max_drain/2,col=hand53_ream_01_sel_seg$max_drain,pch=19)
      
      #plot decli
      plot(hand53_ream_m0_sel_seg$hand,fern_mds$MDS1[crit],cex=hand53_ream_m0_sel_seg$decli/2,pch=19)
      plot(hand53_ream_01_sel_seg$hand,fern_mds$MDS1[crit],cex=hand53_ream_01_sel_seg$decli/2,pch=19)
      

      #plot acum_c
      plot(hand53_ream_m0_sel_seg$hand,fern_mds$MDS1[crit],col=log(hand53_ream_m0_sel_seg$area_c),pch=19,cex=2)
      plot(hand53_ream_01_sel_seg$hand,fern_mds$MDS1[crit],cex=log(hand53_ream_01_sel_seg$area_c)/10,pch=19)
      plot(hand53_ream_01_sel_seg$area_c,fern_mds$MDS1[crit])
      
      hist(hand53_ream_01_sel_seg$area_c/10^4)
      
#multivariate regression tree
data <- cbind(fern_mds[,c("MDS1","MDS2")],hand53_ream_m0_sel[,c("hand","area_c","decli","drain","max_drain")])
data <- decostand(data,method="standardize",MARGIN = 2)
tree <- mvpart(MDS1~hand+max_drain+area_c+decli,data=data)
  summary(tree)
  
  #########################################
  #Moisture index
  MIferns <- read.csv("MIferns_persubunit_raw.csv")
  head(MIferns)
  
  k=752
  MIferns_sel <- subset(MIferns,MIferns$TrNumber==k)
  MI <- rowSums(MIferns_sel[,-c(1,2)])
  topo_sel2 <- subset(topo_sel,topo_sel$trnumber==k)
  
  plot(topo_sel$topo,type="l")
  par(new=T)
 plot(dist,MI)
  dist
  
  #########################################
  
  #load the topography data
  topography <- read.csv("topography.csv", header = TRUE, sep = ",", dec = ".", na.strings = NA)
  head(topography)
  
  #k=786
  #interpolate the topography
  listtopointerpolation <- list()
  for (k in unique(trvector)){
    topo_sel <- subset(topography, topography$tr==k)
    topoapprox <- approx(topo_sel$sub*25,topo_sel$topo,dist)
    plot(topo_sel$topo,type="l",cex=3)
    par(new=T)
    plot(dist)
    
    matrixtopo <- matrix(NA, ncol=3, nrow = length(topoapprox$y))
    matrixtopo[,1] <- topoapprox$y
    matrixtopo[,2] <- rep(k, times=length(topoapprox$y))
    matrixtopo[,3] <- rep(1:20)
    numerator <- which(trvector==k)
    listtopointerpolation[[numerator]] <- matrixtopo
  }
  unlisttopo <- ldply(listtopointerpolation,data.frame)
  
  topo_sel <- unlisttopo[-exx,]
  colnames(topo_sel) <- c("topo","trnumber","sub")
  
  
  
  #solo e hand
  solo_hand <- solo[solo$TrNumber %in% trvector,]
  solo_hand <- solo_hand[-exx,]
  solo_hand$Truni <- as.numeric(paste(solo_hand$TrNumber,solo_hand$subunit,sep = ""))
  hand53_ream_m0_sel$Truni <- as.numeric(paste(hand53_ream_m0_sel$tr,hand53_ream_m0_sel$sub,sep=""))
  fern_mds$Truni <- as.numeric(paste(fern_mds$TrNumber,fern_mds$subunit,sep=""))
  topo_sel$Truni <- as.numeric(paste(topo_sel$trnumber,topo_sel$sub,sep=""))
  
  for (i in solo_hand$Truni){
    solo_hand[solo_hand$Truni==i, "hand" ] <-  hand53_ream_m0_sel[hand53_ream_m0_sel$Truni==i, "hand"]
    solo_hand[solo_hand$Truni==i, "MDS1" ] <-  fern_mds[fern_mds$Truni==i, "MDS1"]
    solo_hand[solo_hand$Truni==i, "MDS2" ] <-  fern_mds[fern_mds$Truni==i, "MDS2"]
    solo_hand[solo_hand$Truni==i, "topo" ] <-  topo_sel[topo_sel$Truni==i, "topo"]
  }

  plot(log(solo_hand$sum_of_basis),-1*solo_hand$MDS1)
  plot(solo_hand$hand,solo_hand$MDS1)
  x <-log(solo_hand$sum_of_basis)
  x <-solo_hand$hand
    y<- solo_hand$MDS1
    
    plot(solo_hand$hand,log(solo_hand$sum_of_basis))
    text(solo_hand$hand,log(solo_hand$sum_of_basis),solo_hand$TrNumber)
    
    regatbc <- lm(solo_hand$MDS1~log(solo_hand$sum_of_basis+solo_hand$hand));summary(regatbc)
  
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
  
##########################################################################################
  #mixed models
  
  mm <- lmer(MDS1 ~ hand + sum_of_basis + (1|TrNumber),solo_hand, REML=FALSE)
  summary(mm)
  
  mmprof <- profile(mm)
  xyplot(mmprof)
  splom(mmprof)
  
  dotplot(ranef(mm, postVar = TRUE))
  qqmath(ranef(mm, postVar=TRUE))
  
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
  
  
  
  reghand <- lm(hand53_ream_m0_sel~scores(mds.ab)[,1]+scores(mds.ab)[,2]+scores(mds.ab)[,3]+scores(mds.ab)[,4]+scores(mds.ab)[,5])
    summary(reghand)  
  
  #regressões
  
  
  #hand e solos
  reg1 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~hand53_ream_m0_sel+log(solo_sel$sum_of_base));summary(reg1)
  reg2 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~hand53_ream_01_sel+solo_sel$sum_of_base);summary(reg2)
  
  #interaction hand e solo    
  reg3 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~hand53_ream_m0_sel*log(solo_sel$sum_of_base));summary(reg3)
  reg4 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~hand53_ream_01_sel*log(solo_sel$sum_of_base));summary(reg4)
  
  #hand only  
  reg5 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~hand53_ream_m0_sel);summary(reg5)
  reg6 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~hand53_ream_m0_sel);summary(reg6)
  
  #topo only  
  reg7 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo);summary(reg7)
  #topo e solo
  reg8 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo+log(solo_sel$sum_of_base));summary(reg8)
  #topo e solo interaction
  reg9 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo*log(solo_sel$sum_of_base));summary(reg9)
  
  #topo hand e solo
  reg10 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo+log(solo_sel$sum_of_base)+hand53_ream_m0_sel);summary(reg10)
  reg11 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo+log(solo_sel$sum_of_base)+hand53_ream_01_sel);summary(reg11)
  #interaction topo hand e solo
  reg12 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo*log(solo_sel$sum_of_base)*hand53_ream_m0_sel);summary(reg12)
  reg13 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo*log(solo_sel$sum_of_base)*hand53_ream_01_sel);summary(reg13)
  
  #topo e hand (exclui solo)      
  reg14 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo+hand53_ream_m0_sel);summary(reg14)
  reg15 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo+hand53_ream_01_sel);summary(reg15)
  #topo e hand interaction
  reg16 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo*hand53_ream_m0_sel);summary(reg16)
  reg17 <- lm(scores(mds.ab)[,1]+scores(mds.ab)[,2]~topo_sel$topo*hand53_ream_01_sel);summary(reg17)
  
  #criar uma tabela de output
  output_reg <- matrix(NA, ncol = 5, nrow = 17)
  colnames(output_reg) <- c("adjR2","AIC","VIFmax","tr","model")
  rownames(output_reg) <- c("hand_solos_m0","hand_solos_01","int_hand_solos_m0","int_hand_solos_01",
                            "hand_m0","hand_01","topo","topo_solo","int_topo_solo","topo_solo_hand_m0",
                            "topo_solo_hand_01","int_topo_hand_solo_m0","int_topo_hand_solo_01",
                            "topo_hand_m0", "topo_hand_01","int_topo_hand_m0","int_topo_hand_01") 
  for (i in 1:17){
    a <- paste("reg",i,sep="")
    output_reg[i,c(1:2)] <- as.numeric(c(round(summary(get(a))$adj.r.squared,2),round(AIC(get(a)),1)))
  }
  for (i in 1:4){
    a <- paste("reg",i,sep="")
    output_reg[i,3] <- round(max(vif(get(a))),2)
  }
  for(i in 8:17){
    a <- paste("reg",i,sep="")
    output_reg[i,3] <- round(max(vif(get(a))),2)
  }
  for(i in 1:17){
    output_reg[i,4] <- k
  }
  output_reg[,5] <- as.character(rownames(output_reg))
  
  #mixed models
  
  geoID_sel <- geoID[geoID$transect %in% trvector,] 
  str(geoID_sel)
  reltopodiff_sel <- reltopodiff[names(reltopodiff)%in%  trvector]
  reltopodiff_sel <- as.character(reltopodiff_sel)
  
  data <- data.frame(scores(mds.ab)[,1],scores(mds.ab)[,2],topo_sel$topo,log(solo_sel$sum_of_base),hand53_ream_m0_sel,hand53_ream_01_sel,geoID_sel$surface,reltopodiff_sel)
  
  colnames(data) <- c("mds1","mds2","topo","log_sb","handm0","hand01","surface","reldiff")
  mm1 <- lmer(mds1+mds2~topo+log_sb+handm0+topo*log_sb*handm0 + (1|surface)+(1|reldiff),data)
  mm2 <- lmer(mds1+mds2~topo+log_sb+handm0+topo*log_sb*handm0 + (1|surface),data)
  mm3 <- lmer(mds1+mds2~topo+log_sb+handm0+topo*log_sb*handm0 + (1|reldiff),data)
  anova(mm2,mm3)
###################################################################################################  
#transfor the list in a data.frame

combineoutputreg <- ldply(listresultsregressions,data.frame)

#fu** the factors
write.csv(combineoutputreg,"lista_teste.csv",row.names = FALSE)
combineoutputreg  <- read.csv("lista_teste.csv",stringsAsFactors = FALSE)


#plot the results
g1 <- ggplot(combineoutputreg,aes(model,adjR2))+ 
  geom_boxplot()
g1


#####################################################################################################
