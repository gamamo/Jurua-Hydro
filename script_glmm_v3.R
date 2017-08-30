
# GLMM models PA

# first eliminate all species with frequency less than 20 occurrences
      
        freqlist <- list()
        for (i in 1:length(speciesall_tomodel_pa)){
          freqlist[[i]] <- data.frame(f = length(which(speciesall_tomodel_pa[,i]!=0)),
                                    sp = i)
        }
        freqdf <- ldply(freqlist,data.frame)
        freqdf <- freqdf[which(freqdf$f>20),]
  
        speciesall_tomodel_pa <- speciesall_tomodel_pa[,freqdf$sp]
        
        speciesall_tomodel_pa <- cbind(speciesall_tomodel_pa,speciesall[,c(1:2)])
        

# join the species and the envi for the modelling and scale variables
        
   
    datamodelsPA <- join(speciesall_tomodel_pa, envi,by = c("TrNumber","subunit"), type = "inner", match = "all" )  
    
# species occorrence in geological surfaces
    colnames(datamodelsPA)
    
    melt1 <- melt(datamodelsPA[,c(1:190,202)])
    melt1 <- melt1[which(melt1$value!=0),]
    spgeo <- dcast(melt1,variable~surface,value.var = "value")
      spgeo[which(spgeo$`Içá formation`!=0),"Içá formation"] <-1
      spgeo[which(spgeo$`Rivers Terraces`!=0),"Rivers Terraces"] <-1
      spgeo[which(spgeo$`Solimões formation`!=0),"Solimões formation"] <-1
      colnames(spgeo)[1] <- "Species"
    
#scale the variables and log trasform them
      
    datamodelsPA$soloslinearint <- log10(datamodelsPA$soloslinearint)
    datamodelsPA$hand <- log10(1+datamodelsPA$hand)
    colnames(datamodelsPA)
    
    toscale <- c("hand","decli","drain","soloslinearint", "solostopoint","srtm_alt")
    datamodelsPA[,toscale] <- scale(datamodelsPA[,toscale])
    datamodelsPA$TrNumber <- as.factor(datamodelsPA$TrNumber)
    datamodelsPA$surface  <- as.factor(datamodelsPA$surface)
    
    
# model 
    colnames(datamodelsPA)
    length(unique(datamodelsPA$TrNumber))
    length(datamodelsPA$subunit)
    seq=seq(1,190)


    listoutputmodelsPA <- list()
    for(i in seq){
      
        m1 <- glmer(datamodelsPA[,i]~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
        m2 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + poly(drain,2) + poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
        m3 <- glmer(datamodelsPA[,i]~poly(hand,2) + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial(link = "logit"),data=datamodelsPA)
        m4 <- glmer(datamodelsPA[,i]~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
        #m5 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + drain  + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
        m6 <- glmer(datamodelsPA[,i]~hand + decli + drain + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)
  

      listaics <- list(m1,m2,m3,m4,m6)
      aictable <- AIC(m1,m2,m3,m4,m6)

      listoutputmodelsPA[[i]]   <- listaics[[which.min(aictable$AIC)]]
    }
    
    listoutputGAM <- list()
    for(k in seq){
      listoutputGAM[[k]] <- gam(datamodelsPA[,k]~s(hand) + s(decli) + drain + s(soloslinearint) ,family = binomial,data=datamodelsPA)
    }

###################################################################    

visreg(a,scale = "response")  
i="Trich.pinn"

m1 <- glmer(Metax.park~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA,glmerControl(optimizer="nloptwrap"))
m2 <- glmer(Metax.park~poly(hand,2) + poly(decli,2) + poly(drain,2) + poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
m3 <- glmer(Metax.park~poly(hand,2) + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial(link = "logit"),data=datamodelsPA)
m4 <- glmer(Metax.park~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
m5 <- glmer(Metax.park~poly(hand,2) + poly(decli,2) + drain  + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
m6 <- glmer(Metax.park~hand + decli + drain + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)

update(m1)


AIC(m1,m2,m3,m4,m5,m6)
tt <- getME(m1,"theta")
ll <- getME(m1,"lower")
min(tt[ll==0])

ss <- getME(m1,c("theta","fixef"))
m1.2 <- update(m1,start=ss,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e4)))

    visreg(m4,scale="response")
 summary(m1)
    sjp.glmer(a)
    sjp.glmer(a2,type = "fe")
    sjp.glmer(m2,type = "eff")
    sjp.glmer(a6)
    
    sjp.glmer(a2,type = "ri.slope")
    sjp.glmer(a2,type = "fe.slope")
    sjp.glmer(a,type = "re.qq")
    
  
  
    visreg(a,"soloslinearint",scale="response",rug=2,alpha=0.001)
    visreg(a,"soloslinearint",rug=2)

##########################################################################################

#obtain some coefficients    
    
    listoutputmodelsPA_table <- list()
    
    for(i in seq){
      listoutputmodelsPA_table[[i]] <- data.frame(Variable = rownames(summary(listoutputmodelsPA[[i]])$coef),
                                               Coefficient = summary(listoutputmodelsPA[[i]])$coef[, 1],
                                               SE = summary(listoutputmodelsPA[[i]])$coef[, 2],
                                               p = round(summary(listoutputmodelsPA[[i]])$coef[, 4],3),
                                               modelName = colnames(datamodelsPA)[i])
    }
    
    
   
    
    toplotPAmodels <- ldply(listoutputmodelsPA_table,data.frame)
      toplotPAmodels[,5]
      toplotPAmodels <- toplotPAmodels[-which(toplotPAmodels$modelName=="indet"),] # to discard the indets in the zingib file
    speciesgroups <- c(rep("Ferns",334),rep("Zingiberales",382),rep("Palms",285),rep("Melastomataceae",217))
    toplotPAmodels <-cbind(toplotPAmodels,speciesgroups)
    toplotPAmodels[order(toplotPAmodels$Coefficient),]
    

    #how many species of each group were modelled?
    
    f <- subset(toplotPAmodels,toplotPAmodels$speciesgroups=="Ferns")
      length(unique(f$modelName))
    z <- subset(toplotPAmodels,toplotPAmodels$speciesgroups=="Zingiberales")
      length(unique(z$modelName))
    p <- subset(toplotPAmodels,toplotPAmodels$speciesgroups=="Palms")
    length(unique(p$modelName))
    m <- subset(toplotPAmodels,toplotPAmodels$speciesgroups=="Melastomataceae")
    length(unique(m$modelName))
  
    #select only the values significant
    
    datagraph_sub <- subset(toplotPAmodels, toplotPAmodels$Variable!="(Intercept)")
    datagraph_sub <- datagraph_sub[which(datagraph_sub$p<0.05),]
    datagraph_sub$p <- round(datagraph_sub$p,3)
    
#how many species have soil significant and how many have hidro variables?
    if(F){
    soilsig <- subset(datagraph_sub,datagraph_sub$Variable=="soloslinearint")
    soilsig$effect <- NA
    soilsig[soilsig$Coefficient>0,"effect"] <- "positive"
    soilsig[soilsig$Coefficient<0,"effect"] <- "negative"
    soilsigdf <- soilsig[order(soilsig$speciesgroups,soilsig$effect,soilsig$modelName),c("modelName","speciesgroups","effect")]
    
    handsig <- subset(datagraph_sub,datagraph_sub$Variable=="hand")
    handsig$effect <- NA
    handsig[handsig$Coefficient>0,"effect"] <- "positive"
    handsig[handsig$Coefficient<0,"effect"] <- "negative"
    handsigdf <- handsig[order(handsig$speciesgroups,handsig$effect,handsig$modelName),c("modelName","speciesgroups","effect")]
    
    declisig <- subset(datagraph_sub,datagraph_sub$Variable=="decli")
    declisig$effect <- NA
    declisig[declisig$Coefficient>0,"effect"] <- "positive"
    declisig[declisig$Coefficient<0,"effect"] <- "negative"
    declisigdf <- declisig[order(declisig$speciesgroups,declisig$effect,declisig$modelName),c("modelName","speciesgroups","effect")]
    
    drainsig <- subset(datagraph_sub,datagraph_sub$Variable=="drain")
    drainsig$effect <- NA
    drainsig[drainsig$Coefficient>0,"effect"] <- "positive"
    drainsig[drainsig$Coefficient<0,"effect"] <- "negative"
    drainsigdf <- drainsig[order(drainsig$speciesgroups,drainsig$effect,drainsig$modelName),c("modelName","speciesgroups","effect")]
    
    }
    ###
#########################################################################################
#### Venn diagrams
    
    
    fernsVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Ferns")
    levels(fernsVenn$Variable)
    fernsVenn$Variable <- revalue(fernsVenn$Variable, c("decli"="Slope","drain"="Drainage","hand"="HAND",
                                                        "poly(soloslinearint, 2)1"="Cation","poly(soloslinearint, 2)2"="Cation",
                                                        "poly(hand, 2)1"="HAND","poly(hand, 2)2"="HAND","poly(decli, 2)1"="Slope",
                                                        "poly(decli, 2)2"="Slope","poly(drain, 2)1"="Drainage","poly(drain, 2)2"="Drainage",
                                                        "soloslinearint"="Cation"))
    
    fernsVennWF <- dcast(fernsVenn,modelName~Variable,fun.aggregate = length)
    colnames(fernsVennWF)[1] <-"Species"
    fernsTiles <- melt(fernsVennWF)
    fernsTiles[fernsTiles$value==2,"value"] <- 1
    
        fernsTilesplot <- ggplot(fernsTiles, aes(variable,Species))+
                                   geom_tile(aes(fill=value), colour = "white")+
                                   scale_fill_gradient(low = "white",high = "steelblue")+
                                   theme(legend.position="none")+
                                   ggtitle("Ferns")
        fernsTilesplot
    
    fernsresuGEO   <- spgeo[spgeo$Species %in% fernsVennWF$Species,]
    fernsresuGEOmelt   <- melt(fernsresuGEO)
    fernsresuGEOmelt$variable <- revalue(fernsresuGEOmelt$variable, c("Içá formation"="IF","Rivers Terraces"="RT","Solimões formation"="SF"))

    
    fernsTilesplotgeo <- ggplot(fernsresuGEOmelt, aes(variable,Species))+
      geom_tile(aes(fill=value), colour = "white")+
      scale_fill_gradient(low = "white",high = "darkolivegreen4")+
      theme(legend.position="none")+
      ggtitle("")+
      ylab("")+
      xlab("")+
      theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())
    fernsTilesplotgeo
    
    
     fga <- ggarrange(fernsTilesplot,fernsTilesplotgeo,widths = c(4, 1))
    
    
   ####### 
    zingVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Zingiberales")
    zingVenn$Variable <- revalue(zingVenn$Variable, c("decli"="Slope","drain"="Drainage","hand"="HAND",
                                                        "poly(soloslinearint, 2)1"="Cation","poly(soloslinearint, 2)2"="Cation",
                                                        "poly(hand, 2)1"="HAND","poly(hand, 2)2"="HAND","poly(decli, 2)1"="Slope",
                                                        "poly(decli, 2)2"="Slope","poly(drain, 2)1"="Drainage","poly(drain, 2)2"="Drainage"))
    zingVennWF <- dcast(zingVenn,modelName~Variable,fun.aggregate = length)
    colnames(zingVennWF)[1] <- "Species"
    zingTiles <- melt(zingVennWF)
    zingTiles[zingTiles$value==2,"value"] <- 1
    head(zingTiles)
    
        zingTilesplot <- ggplot(zingTiles, aes(variable,Species))+
                            geom_tile(aes(fill=value), colour = "white")+
                            scale_fill_gradient(low = "white",high = "steelblue")+
                            theme(legend.position="none")+
                            ggtitle("Zingiberales")
        zingTilesplot
          
    
        
        zingresuGEO   <- spgeo[spgeo$Species %in% zingVennWF$Species,]
        zingresuGEOmelt   <- melt(zingresuGEO )
        zingresuGEOmelt$variable <- revalue(zingresuGEOmelt$variable, c("Içá formation"="IF","Rivers Terraces"="RT","Solimões formation"="SF"))
        
        
        zingTilesplotgeo <- ggplot(zingresuGEOmelt, aes(variable,Species))+
          geom_tile(aes(fill=value), colour = "white")+
          scale_fill_gradient(low = "white",high = "darkolivegreen4")+
          theme(legend.position="none")+
          ggtitle("")+
          ylab("")+
          xlab("")+
          theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())
        zingTilesplotgeo
        
        
        zga <- ggarrange(zingTilesplot,zingTilesplotgeo,widths = c(4, 1))
    
   ######### 
    palmsVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Palms")
    palmsVenn$Variable <- revalue(palmsVenn$Variable, c("decli"="Slope","drain"="Drainage","hand"="HAND",
                                                     "poly(soloslinearint, 2)1"="Cation","poly(soloslinearint, 2)2"="Cation",
                                                     "poly(hand, 2)1"="HAND","poly(hand, 2)2"="HAND","poly(decli, 2)1"="Slope",
                                                     "poly(decli, 2)2"="Slope","poly(drain, 2)1"="Drainage","poly(drain, 2)2"="Drainage"))
    palmsVennWF <- dcast(palmsVenn,modelName~Variable,fun.aggregate = length)
    colnames(palmsVennWF)[1] <- "Species"
    palmTiles <- melt(palmsVennWF)
    palmTiles[palmTiles$value==2,"value"] <- 1
    
    
        palmTilesplot <- ggplot(palmTiles, aes(variable,Species))+
                          geom_tile(aes(fill=value), colour = "white")+
                          scale_fill_gradient(low = "white",high = "steelblue")+
                          theme(legend.position="none")+
                          ggtitle("Palms")
        palmTilesplot
    
        palmresuGEO   <- spgeo[spgeo$Species %in% palmsVennWF$Species,]
        palmresuGEOmelt   <- melt(palmresuGEO )
        palmresuGEOmelt$variable <- revalue(palmresuGEOmelt$variable, c("Içá formation"="IF","Rivers Terraces"="RT","Solimões formation"="SF"))
        
        
        palmTilesplotgeo <- ggplot(palmresuGEOmelt, aes(variable,Species))+
          geom_tile(aes(fill=value), colour = "white")+
          scale_fill_gradient(low = "white",high = "darkolivegreen4")+
          theme(legend.position="none")+
          ggtitle("")+
          ylab("")+
          xlab("")+
          theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())
        palmTilesplotgeo
        
        
        pga <- ggarrange(palmTilesplot,palmTilesplotgeo,widths = c(4, 1))
   
  ####### 
    melasVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Melastomataceae")
    melasVenn$Variable <- revalue(melasVenn$Variable, c("decli"="Slope","drain"="Drainage","hand"="HAND",
                                                        "poly(soloslinearint, 2)1"="Cation","poly(soloslinearint, 2)2"="Cation",
                                                        "poly(hand, 2)1"="HAND","poly(hand, 2)2"="HAND","poly(decli, 2)1"="Slope",
                                                        "poly(decli, 2)2"="Slope","poly(drain, 2)1"="Drainage","poly(drain, 2)2"="Drainage"))
    
    melasVennWF <- dcast(melasVenn,modelName~Variable,fun.aggregate = length)
    colnames(melasVennWF)[1] <- "Species"
    melasTiles <- melt(melasVennWF)
    melasTiles[melasTiles$value==2,"value"] <- 1
    
    
          melasTilesplot <- ggplot(melasTiles, aes(variable,Species))+
                            geom_tile(aes(fill=value), colour = "white")+
                            scale_fill_gradient(low = "white",high = "steelblue")+
                            theme(legend.position="none")+
                            ggtitle("Melastomataceae")
          melasTilesplot
    
    
          melasresuGEO   <- spgeo[spgeo$Species %in% melasVennWF$Species,]
          melasresuGEOmelt   <- melt(melasresuGEO )
          melasresuGEOmelt$variable <- revalue(melasresuGEOmelt$variable, c("Içá formation"="IF","Rivers Terraces"="RT","Solimões formation"="SF"))
          
          
          melasTilesplotgeo <- ggplot(melasresuGEOmelt, aes(variable,Species))+
            geom_tile(aes(fill=value), colour = "white")+
            scale_fill_gradient(low = "white",high = "darkolivegreen4")+
            theme(legend.position="none")+
            ggtitle("")+
            ylab("")+
            xlab("")+
            theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())
          melasTilesplotgeo
          
          
          mga <- ggarrange(melasTilesplot,melasTilesplotgeo,widths = c(4, 1))
    
    
    ## tiles plots
    
    tiff(filename = "TilePlots_geo.tiff", height = 12, width = 11,units = "in",res = 180)
    multiplot(fga,zga,pga,mga,cols = 2)
    #multiplot(fernsTilesplot,zingTilesplot,palmTilesplot,melasTilesplot,cols = 2)
    dev.off()
    
    #write the full table to SI
    write.csv(fernsresuGEO, "output_fixedterms_geo_FERNS.csv",row.names = FALSE)
    write.csv(zingresuGEO, "output_fixedterms_geo_ZING.csv",row.names = FALSE)
    write.csv(palmsresuGEO,"output_fixedterms_geo_PALMS.csv",row.names = FALSE )
    write.csv(melasresuGEO,"output_fixedterms_geo_MELAS.csv",row.names = FALSE )
    
    resugeoALL <- rbind(fernsresuGEO,zingresuGEO,palmsresuGEO,melasresuGEO)
    write.csv(resugeoALL, "output_fixedterms_nogeo.csv",row.names = FALSE)
    
    listVenn <- list(fernsVennWF,zingVennWF,palmsVennWF,melasVennWF)
    
    listVennplots <- list()
    for(i in seq(1,4)){
      cross.area=nrow(subset(listVenn[[i]],listVenn[[i]]$HAND==1|listVenn[[i]]$Hand_2==1 &
                               listVenn[[i]]$Cation_1==1|listVenn[[i]]$Cation_2==1))
      area1 = nrow(subset(listVenn[[i]],listVenn[[i]]$HAND==1|listVenn[[i]]$Hand_2==1)) + cross.area
      area2 = nrow(subset(listVenn[[i]],listVenn[[i]]$Cation_1==1|listVenn[[i]]$Cation_2==1)) + cross.area
      
      listVennplots[[i]] <-draw.pairwise.venn(area1, area2, cross.area, 
                                              cat.pos=c(315,15),cex=1,cat.cex=1.5,ext.pos = c(15,15),
                                              print.mode=c("raw","percent"),ind=TRUE)
    }
    
  i=1
   plot(veeneu)
    tiff(filename = "Venn_diagram.tiff", height = 7, width = 7.5,units = "in",res = 180) # set the pdf file
    grid.arrange(gTree(children=listVennplots[[1]]),gTree(children=listVennplots[[2]]),
             gTree(children=listVennplots[[3]]),gTree(children=listVennplots[[4]]),nrow=2,ncol=2)
    dev.off()
i=1

    aa <- melt(listVenn[[i]])
ggplot(aa, aes(variable,species))+
  geom_tile(aes(fill=value), colour = "white")+
  scale_fill_gradient(low = "white",high = "steelblue")

    
 
#########################################################################################    
    
  
  

