
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
      spgeo[which(spgeo$`River terraces`!=0),"River terraces"] <-1
      spgeo[which(spgeo$`Solimões formation`!=0),"Solimões formation"] <-1
      colnames(spgeo)[1] <- "Species"
    
#scale the variables and log trasform them
      
    datamodelsPA$soloslinearint <- log10(datamodelsPA$soloslinearint)
    datamodelsPA$hand <- log10(1+datamodelsPA$hand)
    colnames(datamodelsPA)
    
    toscale <- c("hand","decli","drain","soloslinearint", "solostopoint","srtm_alt")
    datamodelsPA[,toscale] <- scale(datamodelsPA[,toscale],center = TRUE, scale = TRUE)
    datamodelsPA$TrNumber <- as.factor(datamodelsPA$TrNumber)
    datamodelsPA$surface  <- as.factor(datamodelsPA$surface)
    
    
# model 
    colnames(datamodelsPA)
    length(unique(datamodelsPA$TrNumber))
    length(datamodelsPA$subunit)
    seq=seq(1,190)


    listoutputmodelsPA <- list()
    for(i in seq){
      
        m1 <- glmer(datamodelsPA[,i]~hand + decli + drain + poly(soloslinearint,2,raw=TRUE) + (1|TrNumber),
                    family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                             optCtrl=list(maxfun=2e4)))
        m2 <- glmer(datamodelsPA[,i]~poly(hand,2,raw=TRUE) + poly(decli,2,raw=TRUE) + poly(drain,2,raw=TRUE) + poly(soloslinearint,2,raw=TRUE)+(1|TrNumber),
                    family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                             optCtrl=list(maxfun=2e4)))
        m3 <- glmer(datamodelsPA[,i]~poly(hand,2,raw=TRUE) + decli + drain + poly(soloslinearint,2,raw=TRUE) + (1|TrNumber),
                    family = binomial(link = "logit"),data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                                             optCtrl=list(maxfun=2e4)))
        m4 <- glmer(datamodelsPA[,i]~hand + decli + drain + poly(soloslinearint,2,raw=TRUE) + (1|TrNumber),
                    family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                             optCtrl=list(maxfun=2e4)))
        m5 <- glmer(datamodelsPA[,i]~poly(hand,2,raw=TRUE) + poly(decli,2,raw=TRUE) + drain  + poly(soloslinearint,2,raw=TRUE) + (1|TrNumber),
                    family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                             optCtrl=list(maxfun=2e4)))
        m6 <- glmer(datamodelsPA[,i]~hand + decli + drain + soloslinearint + (1|TrNumber),
                    family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                             optCtrl=list(maxfun=2e4)))
  

      listaics <- list(m1,m2,m3,m4,m5,m6)
      aictable <- AIC(m1,m2,m3,m4,m5,m6)

      listoutputmodelsPA[[i]]   <- listaics[[which.min(aictable$AIC)]]
    }
    
    #Run GAMS for all species just in case to compare if needed
    
    listoutputGAM <- list()
    for(k in seq){
      listoutputGAM[[k]] <- gam(datamodelsPA[,k]~s(hand) + s(decli) + drain + s(soloslinearint) ,family = binomial,data=datamodelsPA)
    }
    
    #get predicted values GLMM
    
    listpredictGLMM <- list()
 
    for (i in seq(listoutputmodelsPA)){
      listpredictGLMM[[i]] <-  predict(listoutputmodelsPA[[i]],type="response")
    }
    
    # get maximum predicted values
    listMAXpredi <- list()
    listpredicvalue <- list()
    for (i in seq(listpredictGLMM)){
      listMAXpredi[[i]] <- datamodelsPA[which.max(listpredictGLMM[[i]]),c("TrNumber" ,"subunit","hand","soloslinearint")]
      listpredicvalue[[i]] <- as.vector(sort(listpredictGLMM[[i]],decreasing = TRUE)[1])
    }
    
    

    dfmaxpredic <- ldply(listMAXpredi,data.frame)
    dfvalpredic <- ldply(listpredicvalue,data.frame)
      colnames(dfvalpredic) <- "predict"
    speciesnames <- as.vector(colnames(datamodelsPA)[seq]) #set seq when running the models
    dfmaxpredic <- cbind(speciesnames,dfmaxpredic,dfvalpredic)

###################################################################    
#testes


    i="Adian.peti"
    
    
    m1 <- glmer(Danae.b.g~hand + decli + drain + poly(soloslinearint,2,raw=TRUE) + (1|TrNumber),family = binomial,
                data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                       optCtrl=list(maxfun=2e4)))
    
    m2 <- glmer(Danae.b.g~poly(hand,2,raw=TRUE) + poly(decli,2,raw=TRUE) + poly(drain,2,raw=TRUE) + poly(soloslinearint,2,raw=TRUE)+
                  (1|TrNumber),family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                                        optCtrl=list(maxfun=2e4)))
    m3 <- glmer(Danae.b.g~poly(hand,2,raw=TRUE) + decli + drain + poly(soloslinearint,2,raw=TRUE) + (1|TrNumber),
                family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                         optCtrl=list(maxfun=2e4)))
    m4 <- glmer(Danae.b.g~hand + decli + drain + poly(soloslinearint,2,raw=TRUE) + (1|TrNumber),
                family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e4)))
    m5 <- glmer(Danae.b.g~poly(hand,2,raw=TRUE) + poly(decli,2,raw=TRUE) + drain  + poly(soloslinearint,2,raw=TRUE) + (1|TrNumber),
                family = binomial,data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                                         optCtrl=list(maxfun=2e4)))
    m6 <- glmer(Danae.b.g~hand + decli + drain + soloslinearint + (1|TrNumber),family = binomial,
                data=datamodelsPA,control=glmerControl(optimizer="bobyqa",
                                                       optCtrl=list(maxfun=2e4)))
    
    
    aictable <- AIC(m1,m2,m3,m5,m4,m6)
    aictable 
    visreg(m6,scale = "response")  
    
    Anova(m1)
    summary(m6)
    
    
    a <- glmer(Adian.peti~hand+poly(soloslinearint,2)+(1|TrNumber),family=binomial, data = datamodelsPA)
    summary(a)
    plot(a1)
    qqnorm(a,~ranef(.))
    a1 <- glmmTMB(Adian.peti~hand+poly(soloslinearint,2)+(1|TrNumber),family=nbinom1, data = datamodelsPA)
    summary(a1)

    sjp.glmer(m4)
    sjp.glmer(m4,type = "fe")
    sjp.glmer(m4,type = "eff")
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
    speciesgroups <- c(rep("Pteridophyte",337),rep("Zingiberales",394),rep("Palms",286),rep("Melastomataceae",218))
    toplotPAmodels <-cbind(toplotPAmodels,speciesgroups)
    toplotPAmodels[order(toplotPAmodels$Coefficient),]
    
    
    #how many species of each group were modelled?
    
    f <- subset(toplotPAmodels,toplotPAmodels$speciesgroups=="Pteridophyte")
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
    
    a <- subset(toplotPAmodels, toplotPAmodels$Variable=="poly(soloslinearint, 2, raw = TRUE)1"|toplotPAmodels$Variable=="poly(soloslinearint, 2, raw = TRUE)2")
    a <- subset(toplotPAmodels, toplotPAmodels$Variable=="poly(soloslinearint, 2, raw = TRUE)2")
    
    #creat a colunm to indicate if the effect was positive or negative
    
    datagraph_sub$effect <- NA
    datagraph_sub[datagraph_sub$Coefficient>0,"effect"] <- 1
    datagraph_sub[datagraph_sub$Coefficient<0,"effect"] <- -1

    
#########################################################################################
####  plots
    
#change levels of the variables
    
    datagraph_sub$Variable <- revalue(datagraph_sub$Variable, c("decli"="Slope","drain"="Drainage","hand"="HAND",
                                                                "poly(soloslinearint, 2, raw = TRUE)1"="Cations",
                                                                "poly(soloslinearint, 2, raw = TRUE)2"="Cations2",
                                                                "poly(hand, 2, raw = TRUE)1"="HAND",
                                                                "poly(hand, 2, raw = TRUE)2"="HAND2",
                                                                "poly(decli, 2, raw = TRUE)1"="Slope",
                                                                "poly(decli, 2, raw = TRUE)2"="Slope2",
                                                                "poly(drain, 2, raw = TRUE)1"="Drainage",
                                                                "poly(drain, 2, raw = TRUE)2"="Drainage2",
                                                                "soloslinearint"="Cations"))

    #clean the polinomiais. if the species is significant for both polynomiais, only take the second order one
   
      listcoefclean <- list()
      
      for(i in unique(datagraph_sub$speciesgroups)){
      tiles <- subset(datagraph_sub,datagraph_sub$speciesgroups==i)
      
      cc <- subset(tiles,tiles$Variable=="Cations" | tiles$Variable=="Cations2")
      twoOc <- names(which(table(cc$modelName)==2))
      tosep <- cc[cc$modelName %in% twoOc , ]
      toincl <- cc[!cc$modelName %in% twoOc , ]
      tosep2 <- tosep[which(tosep$Variable=="Cations2"),]
      ccc <- rbind(toincl,tosep2)
      
      hc <- subset(tiles,tiles$Variable=="HAND" | tiles$Variable=="HAND2")
      twoOc <- names(which(table(hc$modelName)==2))
      tosep <- hc[hc$modelName %in% twoOc , ]
      toincl <- hc[!hc$modelName %in% twoOc , ]
      tosep2 <- tosep[which(tosep$Variable=="HAND2"),]
      hhc <- rbind(toincl,tosep2)
      
      sc <- subset(tiles,tiles$Variable=="Slope" | tiles$Variable=="Slope2")
      twoOc <- names(which(table(sc$modelName)==2))
      tosep <- sc[sc$modelName %in% twoOc , ]
      toincl <- sc[!sc$modelName %in% twoOc , ]
      tosep2 <- tosep[which(tosep$Variable=="Slope2"),]
      ssc <- rbind(toincl,tosep2)
      
      dc <- subset(tiles,tiles$Variable=="Drainage" | tiles$Variable=="Drainage2")
      twoOc <- names(which(table(dc$modelName)==2))
      tosep <- dc[dc$modelName %in% twoOc , ]
      toincl <- dc[!dc$modelName %in% twoOc , ]
      tosep2 <- tosep[which(tosep$Variable=="Drainage2"),]
      ddc <- rbind(toincl,tosep2)
      
      listcoefclean[[i]] <- rbind(ccc,hhc,ssc,ddc)
      listcoefclean[[i]]$Variable <- revalue(listcoefclean[[i]]$Variable,c("Cations2"="Cations",
                                            "HAND2"="HAND",
                                            "Slope2"="Slope",
                                            "Drainage2"="Drainage"))
      }
      

#plots
      
      listggplot <- list()
      for (j in names(listcoefclean)){
        
      Tilesplot <- ggplot(listcoefclean[[j]] , aes(Variable,modelName,size= abs(Coefficient)))+
        geom_point(aes(color=effect))+
        scale_colour_gradient2()+
        ggtitle(j)+
        xlab("")+
        ylab("")+
        theme_light()+
        scale_size(range = c(1,10))+
        theme(legend.position="none")+
        theme(axis.text.x = element_text(size = 14))
      
      
      resuGEO   <- spgeo[spgeo$Species %in% listcoefclean[[j]]$modelName,]
      resuGEOmelt   <- melt(resuGEO)
      resuGEOmelt$variable <- revalue(resuGEOmelt$variable, c("Içá formation"="IF","River terraces"="RT","Solimões formation"="SF"))
      
      
      Tilesplotgeo <- ggplot(resuGEOmelt, aes(variable,Species))+
        geom_tile(aes(fill=value), colour = "white")+
        scale_fill_gradient(low = "white",high = "gray30")+
        theme(legend.position="none")+
        ggtitle("")+
        ylab("")+
        xlab("")+
        theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())+
        theme(axis.text.x = element_text(size = 12))

      listggplot[[j]]<- ggarrange(Tilesplot,Tilesplotgeo,widths = c(4, 1))
      
      }

    
    ## tiles plots
    
    tiff(filename = "TilePlots_geo.tiff", height = 12, width = 11,units = "in",res = 96)
    multiplot(listggplot[["Pteridophyte"]],listggplot[["Zingiberales"]],
              listggplot[["Palms"]],listggplot[["Melastomataceae"]],cols = 2)
    dev.off()
    
    #write the full table to SI
    write.csv(fernsresuGEO, "output_fixedterms_geo_FERNS.csv",row.names = FALSE)
    write.csv(zingresuGEO, "output_fixedterms_geo_ZING.csv",row.names = FALSE)
    write.csv(palmsresuGEO,"output_fixedterms_geo_PALMS.csv",row.names = FALSE )
    write.csv(melasresuGEO,"output_fixedterms_geo_MELAS.csv",row.names = FALSE )
    
    resugeoALL <- rbind(fernsresuGEO,zingresuGEO,palmsresuGEO,melasresuGEO)
    write.csv(resugeoALL, "output_fixedterms_nogeo.csv",row.names = FALSE)
    
    
    
    
#########################################################################################
# plot the predicted values
    
    unlistcoefclean <- ldply(listcoefclean,data.frame)
  
    
    datapredictplotWF <- dcast(unlistcoefclean,modelName+speciesgroups~Variable,fun.aggregate = length)
    colnames(datapredictplotWF)[1] <-"Species"
    
    preparePredictdata <- datapredictplotWF[which(datapredictplotWF$HAND>=1 & datapredictplotWF$Cation>=1),]
    
    #select only species where hand and cations were significative
    #this objects were taken formt the beggining of the script
    
    predictplot <- dfmaxpredic[dfmaxpredic$speciesnames %in% preparePredictdata$Species,]
    predictplot <- cbind(predictplot,preparePredictdata$speciesgroups)
    colnames(predictplot)[7] <- "groups"

    colnames(predictplot)

    #redblue <- colorRampPalette(c("red","blue"),bias=.1,space="rgb")
    
    Pplot <- ggplot(predictplot,aes(x=soloslinearint,y=hand,label=speciesnames))+
             geom_point(aes(color=groups))+
             geom_point(aes(size=3,color=groups),show.legend=FALSE)+ 
             #geom_point(aes(size=3),show.legend=FALSE)+
             #geom_point(aes(color=groups),show.legend=TRUE)+
             theme_bw()+
             #scale_color_manual(values = redblue(4))+
             #scale_color_gradient(name="groups",low="blue",high="red")+
             geom_text_repel(size=5)+
             #theme(legend.position="none")+
             theme(text = element_text(size=16))+
             xlab("Soil Cation Concentration")+
             ylab("HAND")+
             theme(legend.position = c(0.9,0.92),legend.title = element_blank(),legend.text = element_text(size = 12)) +
             theme(axis.ticks = element_blank(),axis.text = element_blank())
    Pplot
    
    tiff(filename = "predic_plot.tiff", height = 8, width = 9.5,units = "in",res = 180)
    Pplot
    dev.off()
    
    
    
    
        
#########################################################################################
#not in use    
    
    
    
    
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
    
  
  

