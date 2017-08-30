
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
    seq=seq(1,190)
    

    listoutputmodelsPA <- list()
    for(i in seq){
      
      if(i!=60){
      m1 <- glmer(datamodelsPA[,i]~hand + decli + drain + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)
      m2 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + poly(drain,2) + poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
      m3 <- glmer(datamodelsPA[,i]~poly(hand,2) + decli + drain + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)
      m4 <- glmer(datamodelsPA[,i]~poly(hand,2) + decli + drain + poly(soloslinearint)+(1|TrNumber),family = binomial,data=datamodelsPA)
      m5 <- glmer(datamodelsPA[,i]~hand + decli + drain + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
      m6 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + drain  + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)
      m7 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + drain  + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
      m8 <- gam(datamodelsPA[,i]~s(hand)+s(decli) + drain + s(soloslinearint),family = binomial,data=datamodelsPA)
      } else {
      m1 <- glmer(datamodelsPA[,i]~hand + decli + drain + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)
      m2 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + poly(drain,2) + poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
      m3 <- glmer(datamodelsPA[,i]~poly(hand,2) + decli + drain + soloslinearint + (1|TrNumber),family = binomial,data=datamodelsPA)
      m4 <- glmer(datamodelsPA[,i]~poly(hand,2) + decli + drain + poly(soloslinearint)+(1|TrNumber),family = binomial,data=datamodelsPA)
      m7 <- glmer(datamodelsPA[,i]~poly(hand,2) + poly(decli,2) + drain  + poly(soloslinearint,2) + (1|TrNumber),family = binomial,data=datamodelsPA)
      m8 <- gam(datamodelsPA[,i]~s(hand)+s(decli) + drain + s(soloslinearint),family = binomial,data=datamodelsPA)
      } 
      
      listaics <- list(m1,m2,m3,m4,m5,m6,m7,m8)
      aictable <- AIC(m1,m2,m3,m4,m5,m6,m7,m8)
      
      listoutputmodelsPA[[i]]   <- listaics[[which.min(aictable$AIC)]]
    }
    
    
    visreg(m4,scale="response")
i=60
    summary(m4)
    
    m4 <- glmer(Adian.humi~hand+ decli + drain + poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
    
    a <- glmer(Bolbi.lind~hand+decli+drain+soloslinearint+(1|TrNumber),family = binomial,data=datamodelsPA)
    a1 <- glmer(Bolbi.lind~hand+decli+drain+poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
    a2 <- glmer(Bolbi.lind~poly(hand,2)+decli+drain+poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
    a3 <-  glmer(Bolbi.lind~poly(hand,2)+poly(decli,2)+drain+poly(soloslinearint,2,simple=TRUE)+(1|TrNumber),family = binomial,data=datamodelsPA)
    a4 <-  glmer(Bolbi.lind~poly(hand,2)+poly(decli,2)+poly(drain,2)+poly(soloslinearint,2)+(1|TrNumber),family = binomial,data=datamodelsPA)
    a5 <- glmer(Bolbi.lind~poly(hand,2)+decli+drain+soloslinearint+(1|TrNumber),family = binomial,data=datamodelsPA)
    a6 <- gam(Bolbi.lind~s(hand)+s(decli)+drain+s(soloslinearint),family = binomial,data=datamodelsPA)
    summary(a1)
    AIC(a,a1)
    
    poly(datamodelsPA$soloslinearint,2)
    
    plot(a6)
    plot.gam
    listatest <- list(a,a1,a2,a3,a4,a5)
    
    aictable <- AIC(a,a1,a2,a3,a4,a5,a6)
    listatest[[which.min(aictable$AIC)]]
    
    visreg(a2,scale = "response")
    
###################################################################    
    ## testes
    colnames(datamodelsPA)
    datamodelsPA$solos
    
    a <- glmer(Adian.p.t~hand+decli+drain+soloslinearint+(1|TrNumber),family = binomial,data=datamodelsPA)
    summary(a)
    sort(predict(a,type="response"))
    
    plot(datamodelsPA$soloslinearint,datamodelsPA$hand,
         col=ifelse(rownames(datamodelsPA)=="280"|rownames(datamodelsPA)=="93", "red", "black"))
    
    a <- glmer(Didym.trun~hand+decli+drain+soloslinearint+(1|TrNumber),family = binomial,data=datamodelsPA)
    summary(a)
    sort(predict(a,type="response"))
    
    plot(datamodelsPA$soloslinearint,datamodelsPA$hand,pch=19,
         col=ifelse(rownames(datamodelsPA)=="478"|rownames(datamodelsPA)=="439", "red", "black"))

    
    a <- glmer(Trich.pinn~hand+decli+drain+soloslinearint+(1|TrNumber),family = binomial,data=datamodelsPA)
    summary(a)
    sort(predict(a,type="response"))
    
    plot(datamodelsPA$soloslinearint,datamodelsPA$hand,pch=19,
         col=ifelse(rownames(datamodelsPA)==c("354")|rownames(datamodelsPA)==c("524"), "red", "black"))
    plot(datamodelsPA$drain,datamodelsPA$hand,pch=19,
         col=ifelse(rownames(datamodelsPA)==c("354")|rownames(datamodelsPA)==c("524"), "red", "black"))
    
    
    a <- glmer(Bolbi.lind~hand+decli+drain+soloslinearint+(1|TrNumber),family = binomial,data=datamodelsPA)
    summary(a)
    sort(predict(a,type="response"))
    
    plot(datamodelsPA$soloslinearint,datamodelsPA$hand,pch=19,
         col=ifelse(rownames(datamodelsPA)==c("459")|rownames(datamodelsPA)==c("478"), "red", "black"))
    plot(datamodelsPA$drain,datamodelsPA$hand,pch=19,
         col=ifelse(rownames(datamodelsPA)==c("459")|rownames(datamodelsPA)==c("478"), "red", "black"))
    
    
    
    
 
    sjp.glmer(a2)
    sjp.glmer(a2,type = "fe")
    sjp.glmer(a1,type = "eff")
    sjp.glmer(a6)
    
    sjp.glmer(a2,type = "ri.slope")
    sjp.glmer(a2,type = "fe.slope")
    sjp.glmer(a,type = "re.qq")
    
    sort(predict(a6,type="response"))
    
    var1 = poly(datamodelsPA$hand,2)
    
    hist(poly(datamodelsPA$hand,2)[,2])
    var2 = datamodelsPA$soloslinearint
    
    datateste <- data.frame(var = c("poly(var1,2)","poly(var2,2)"))
    
    model <- glmer(Asple.pear~var1+var2+(1|TrNumber),family = binomial,data=datamodelsPA)
    
    summary(a)
    
    aictable <- AIC(a,a1,a2,a3,a4,a5,a6)
    class(aictable)
    as.symbol(print(rownames(aictable[which.max(aictable$AIC),])))
    
    
    sort(predict(a2,type="response"))
    
    plot(datamodelsPA$soloslinearint,datamodelsPA$hand,pch=19,
         col=ifelse(rownames(datamodelsPA)==c("486")|rownames(datamodelsPA)==c("487"), "red", "black"))
    
    plot(datamodelsPA$drain,datamodelsPA$hand,pch=19,
         col=ifelse(rownames(datamodelsPA)==c("459")|rownames(datamodelsPA)==c("478"), "red", "black"))
    
    visreg(a,scale = "response")
    
    
    
    b <- glmer(Bolbi.nico~hand+decli+drain+soloslinearint+(1|TrNumber)+(1|surface),family = binomial,data=datamodelsPA)
    summary(b)
    
    c <- glmer(Asple.pear~poly(hand,2)+decli+drain+soloslinearint+(1|TrNumber),family = binomial,data=datamodelsPA)
    summary(c)
    Beta(c)
    sort(predict(c,type=c("terms")))
    length(datamodelsPA$Linds.phas)
    plot(c)
    plot(log(1+datamodelsPA$soloslinearint),fitted(c))
    plot(allEffects(c),rug=FALSE)
    visreg(c,)
    
    d <- glmer(Linds.phas~hand+decli+drain+soloslinearint+(surface|TrNumber),family = binomial,data=datamodelsPA)
    e <- glmer(Linds.phas~(surface|TrNumber),family = binomial,data=datamodelsPA)
    
    
    
    f <- glm(Asple.pear~hand+decli+drain+soloslinearint,family = binomial,data=datamodelsPA)
    summary(f)
    sort(predict(f))
    plot(fitted(f),datamodelsPA$soloslinearint)
    plot(f)
    plot(datamodelsPA$soloslinearint,datamodelsPA$hand,cex)
    
    
    g <- gam(Linds.phas~s(hand)+s(soloslinearint),family = binomial,data=datamodelsPA)
    summary(g)
    
    h <- gam(Linds.phas~s(log(1+soloslinearint)),family = binomial,data=datamodelsPA)
    summary(h)
    plot(gam)
    max(log(1+datamodelsPA$soloslinearint))
    
    MyData <- data.frame(soloslinearint = seq(from=0.10, to =2.96, by=0.5))
    Pred <- predict(h,type="response",newdata=MyData)
    
        plot(log(1+datamodelsPA$soloslinearint),datamodelsPA$Linds.phas)
        lines(MyData$soloslinearint,log(1+Pred))

    
    m <- gamlss(Adian.p.t~hand,data=datamodelsPA,family = BI)
    str(g)
    anova(a,g,f)
    AIC(g,f,a)
    
    mgcv::plot.gam(g,pages=1,,residuals=TRUE,all.terms = TRUE,shade=TRUE,shade.col=2)
    plot(g,all.terms=TRUE,se=TRUE)
    plot(f)
    summary(f)
    anova(g)
    gam.check(g)
      
    Ef2 <- resid(f)
    F2  <- fitted(f) 
    plot(F2,Ef2)
    boxplot(Ef2~TrNumber, data=datamodelsPA)
    boxplot(Ef2~surface, data=datamodelsPA)
    abline(0,0)
    plot(datamodelsPA$soloslinearint,F2)
    
   
    
    summary(c)
    Anova(c) 
    
    
    
    AIC(a,b,c,d,e)
    
    i=4
    
    i="lepitenu"
    i="cala.sp32"
    i="Adian.pulv"
    i="Metax.rost"
    
    plot(log(datamodelsPA$soloslinearint),datamodelsPA[,"Linds.phas"])
    aa <- cbind(datamodelsPA$Bolbi.nico,datamodelsPA$sum_of_basis,datamodelsPA$surface,datamodelsPA$soloslinearint)
    aa[order(datamodelsPA$Bolbi.nico,decreasing = TRUE),]
    
    a <- glmer(Bolbi.nico~hand+drain+decli+soloslinearint+(1|surface/TrNumber),family = binomial,data=datamodelsPA)
    b <- glmer(Adian.term~hand+decli+drain+soloslinearint+(1|TrNumber)+(1|surface),family = binomial,data=datamodelsPA)
    c <- glmer(Adian.term~hand+decli+drain+soloslinearint+(1+surface|TrNumber),family = binomial,data=datamodelsPA)
    
    AIC(b,c)
    coefplot(list(a,b))
    coef(a)

    coef(listoutputmodelsPA[[2]])
    a <- listoutputmodelsPA[[2]]
    ranef(listoutputmodelsPA[[2]])
    ranef(a)
    coef(b)$surface
    
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
    
str(a)

GAM_status <- function () {
  if (all(c( "mgcv") %in% .packages())) print("Not OK")
  else print("OK")
}
GAM_status()
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
    speciesgroups <- c(rep("Ferns",255),rep("Zingiberales",300),rep("Palms",220),rep("Melastomataceae",175))
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
    
    
    ###
#########################################################################################
#### Venn diagrams
    
    
    fernsVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Ferns")
    fernsVennWF <- dcast(fernsVenn,modelName~Variable,fun.aggregate = length)
    colnames(fernsVennWF) <- c("species","Slope","Drainage","HAND","Cation")
    
    fernsresuGEO <- join(fernsVennWF, spgeo ,by = c("species"), type = "left", match = "all")
    
    zingVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Zingiberales")
    zingVennWF <- dcast(zingVenn,modelName~Variable,fun.aggregate = length)
    colnames(zingVennWF) <- c("species","Slope","Drainage","HAND","Cation")
    
    zingresuGEO <- join(zingVennWF, spgeo ,by = c("species"), type = "left", match = "all")
    
    palmsVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Palms")
    palmsVennWF <- dcast(palmsVenn,modelName~Variable,fun.aggregate = length)
    colnames(palmsVennWF) <- c("species","Slope","Drainage","HAND","Cation")
    
    palmsresuGEO <- join(palmsVennWF, spgeo ,by = c("species"), type = "left", match = "all")
    
    melasVenn <- subset(datagraph_sub,datagraph_sub$speciesgroups=="Melastomataceae")
    melasVennWF <- dcast(melasVenn,modelName~Variable,fun.aggregate = length)
    colnames(melasVennWF) <- c("species","Slope","HAND","Cation")
    melasVennWF$Drainage <- rep(0, length(melasVennWF$species))
    
    
    melasresuGEO <- join(melasVennWF, spgeo ,by = c("species"), type = "left", match = "all")
    
    #write the full table to SI
    resugeoALL <- rbind(fernsresuGEO,zingresuGEO,palmsresuGEO,melasresuGEO)
    write.csv(resugeoALL, "output_fixedterms_nogeo.csv",row.names = FALSE)
    
    listVenn <- list(fernsVennWF,zingVennWF,palmsVennWF,melasVennWF)
    
    listVennplots <- list()
    for(i in seq(1,4)){
      area1 = nrow(subset(listVenn[[i]],listVenn[[i]]$HAND==1))
      area2 = nrow(subset(listVenn[[i]],listVenn[[i]]$Cation==1))
      cross.area=nrow(subset(listVenn[[i]],listVenn[[i]]$HAND==1 &  listVenn[[i]]$Cation==1))
      listVennplots[[i]] <-draw.pairwise.venn(area1, area2, cross.area,category = c("HAND","Cations"), 
                                              cat.pos=c(315,15),cex=1,cat.cex=1.5,ext.pos = c(15,15),print.mode=c("raw","percent"),ind=FALSE)
    }
    
  
   
    tiff(filename = "Venn_diagram.tiff", height = 7, width = 7.5,units = "in",res = 180) # set the pdf file
    grid.arrange(gTree(children=listVennplots[[1]]),gTree(children=listVennplots[[2]]),
             gTree(children=listVennplots[[3]]),gTree(children=listVennplots[[4]]),nrow=2,ncol=2)
    dev.off()

    
    
 
#########################################################################################    
    
  
  

