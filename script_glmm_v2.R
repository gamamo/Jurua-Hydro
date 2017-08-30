
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
      colnames(spgeo)[1] <- "species"
    
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
    seq=seq(1,190)[-c(38,77)]

    

    listoutputmodelsPA <- list()
    for(i in seq){
      
      
        m1 <- glmer(datamodelsPA[,i]~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
        m2 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + poly(drain,2) + poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
        m3 <- glmer(datamodelsPA[,i]~poly(hand,2) + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial(link = "logit"),data=datamodelsPA)
        m4 <- glmer(datamodelsPA[,i]~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
        m5 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + drain  + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
        m6 <- glmer(datamodelsPA[,i]~hand + decli + drain + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)
  

      listaics <- list(m1,m2,m3,m4,m5)
      aictable <- AIC(m1,m2,m3,m4,m5)

      listoutputmodelsPA[[i]]   <- listaics[[which.min(aictable$AIC)]]
    }
    
    listoutputGAM <- list()
    for(i in seq){
      listoutputGAM[[i]] <- gam(datamodelsPA[,i]~s(hand) + s(decli) + drain + s(soloslinearint) ,family = binomial,data=datamodelsPA)
    }

###################################################################    
a <- gam(Bolbi.nico~s(hand) + s(decli) + drain + s(soloslinearint) ,family = binomial,data=datamodelsPA)
summary(a)

m1 <- glmer(Asple.pear~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
m2 <- glmer(Asple.pear~poly(hand,3) + poly(decli,2) + poly(drain,2) + poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
m3 <- glmer(Asple.pear~poly(hand,2) + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial(link = "logit"),data=datamodelsPA)
m4 <- glmer(Asple.pear~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
m5 <- glmer(Asple.pear~poly(hand,2) + poly(decli,2) + drain  + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
m6 <- glmer(Asple.pear~hand + decli + drain + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)

visreg(m2,scale = "response")    
visreg(a,scale = "response")  

    i="heli.acum"
 summary(m2)
    sjp.glmer(a)
    sjp.glmer(a2,type = "fe")
    sjp.glmer(m2,type = "eff")
    sjp.glmer(a6)
    
    sjp.glmer(a2,type = "ri.slope")
    sjp.glmer(a2,type = "fe.slope")
    sjp.glmer(a,type = "re.qq")
    
    
    visreg(a,scale = "response")
    ?convergence
    
  
  
    summary(a)
    sjp.glmer(a,type = "fe")
    sjp.glmer(a,type = "eff")
    sjp.glmer(g)

    sjp.glmer(a,type = "ri.slope")
    sjp.glmer(a,type = "fe.slope")
    sjp.glmer(a,type = "re.qq")
    
    Ef2 <- resid(a)
    F2  <- fitted(a) 
    plot(F2,Ef2)
    boxplot(Ef2~TrNumber, data=datamodelsPA)
    boxplot(coef(a)$surface$~surface, data=datamodelsPA)
    abline(0,0)
    
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
    
    
    toplotPAmodels <- toplotPAmodels[-which(toplotPAmodels$modelName=="indet"),]
    
    toplotPAmodels <- ldply(listoutputmodelsPA_table,data.frame)
      toplotPAmodels[,5]
    speciesgroups <- c(rep("Ferns",341),rep("Zingiberales",403),rep("Palms",295),rep("Melastomataceae",231))
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
    fernsVennWF <- dcast(fernsVenn,modelName~Variable,fun.aggregate = length)
    colnames(fernsVennWF) <- c("species","Slope","Drainage","HAND","Cation_1","Cation_2","HAND","Hand_2","Slope", "Slope_2","Drainage_2")
    fernsTiles <- melt(fernsVennWF)
    fernsTiles$variable <- revalue(fernsTiles$variable, c("Hand_2"="HAND","Cation_1"="Cation","Cation_2"="Cation",
                                                          "Slope_2"="Slope","Drainage_2"="Drainage"))
        
        fernsTilesplot <- ggplot(fernsTiles, aes(variable,species))+
                                   geom_tile(aes(fill=value), colour = "white")+
                                   scale_fill_gradient(low = "white",high = "steelblue")+
                                   theme(legend.position="none")+
                                   ggtitle("Ferns")
        fernsTilesplot
        
    
    fernsresuGEO <- join(fernsVennWF, spgeo ,by = c("species"), type = "left", match = "all")
    
    zingVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Zingiberales")
    zingVennWF <- dcast(zingVenn,modelName~Variable,fun.aggregate = length)
    colnames(zingVennWF) <- c("species","Slope","Drainage","HAND","Cation_1","Cation_2","HAND","Hand_2","Slope", "Slope_2","Drainage","Drainage_2")
    zingTiles <- melt(zingVennWF)
    zingTiles$variable <- revalue(zingTiles$variable, c("Hand_2"="HAND","Cation_1"="Cation","Cation_2"="Cation",
                                                        "Slope_2"="Slope","Drainage_2"="Drainage"))
    
        zingTilesplot <- ggplot(ZingTiles, aes(variable,species))+
                            geom_tile(aes(fill=value), colour = "white")+
                            scale_fill_gradient(low = "white",high = "steelblue")+
                            theme(legend.position="none")+
                            ggtitle("Zingiberales")
        zingTilesplot
          
    
    zingresuGEO <- join(zingVennWF, spgeo ,by = c("species"), type = "left", match = "all")
    
    palmsVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Palms")
    palmsVennWF <- dcast(palmsVenn,modelName~Variable,fun.aggregate = length)
    colnames(palmsVennWF) <- c("species","Slope","Drainage","HAND","Cation_1","Cation_2","HAND","Hand_2","Slope", "Slope_2","Drainage_2")
    palmTiles <- melt(palmsVennWF)
    palmTiles$variable <- revalue(palmTiles$variable, c("Hand_2"="HAND","Cation_1"="Cation","Cation_2"="Cation",
                                                        "Slope_2"="Slope","Drainage_2"="Drainage"))
    
        palmTilesplot <- ggplot(palmTiles, aes(variable,species))+
                          geom_tile(aes(fill=value), colour = "white")+
                          scale_fill_gradient(low = "white",high = "steelblue")+
                          theme(legend.position="none")+
                          ggtitle("Palms")
        palmTilesplot
    
    palmsresuGEO <- join(palmsVennWF, spgeo ,by = c("species"), type = "left", match = "all")
    
    melasVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Melastomataceae")
    melasVennWF <- dcast(melasVenn,modelName~Variable,fun.aggregate = length)
    colnames(melasVennWF) <- c("species","Slope","Drainage","HAND","Cation_1","Cation_2","HAND","Hand_2", "Slope_2","Drainage_2")
    melasTiles <- melt(melasVennWF)
    melasTiles$variable <- revalue(melasTiles$variable, c("Hand_2"="HAND","Cation_1"="Cation","Cation_2"="Cation",
                                                        "Slope_2"="Slope","Drainage_2"="Drainage"))
    
          melasTilesplot <- ggplot(melasTiles, aes(variable,species))+
                            geom_tile(aes(fill=value), colour = "white")+
                            scale_fill_gradient(low = "white",high = "steelblue")+
                            theme(legend.position="none")+
                            ggtitle("Melastomataceae")
          melasTilesplot
    
    
    melasresuGEO <- join(melasVennWF, spgeo ,by = c("species"), type = "left", match = "all")
    
    
    ## tiles plots
    
    tiff(filename = "TilePlots.tiff", height = 12, width = 8.5,units = "in",res = 180)
    multiplot(fernsTilesplot,zingTilesplot,palmTilesplot,melasTilesplot,cols = 2)
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
    
  
  

