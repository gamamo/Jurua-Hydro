# script to prepare a table that contains the mean optima
# of the species per subunit
# tun the script when all the objects from the scripts:
# "jurua_topoHydro" and "landscape_scale_analysis" have been run


getwd()

dir()
optimaAll <- read.csv("otimosALL_gm.csv",sep=";",stringsAsFactors = FALSE)

  head(optimaAll)
  length(optimaAll[,1])
  
 
  optimaAll$spcode_1 <- substring(optimaAll$spnames,1,5)
  optimaAll$spcode_2 <- substring(optimaAll$epiteto,1,4)
  optimaAll$spcode   <- paste(optimaAll$spcode_1,".",optimaAll$spcode_2,sep = "")
  
  colnames(fern25)
  opt <- optimaAll[optimaAll$spcode %in% colnames(fern25),]
  length(opt[,1])
  head(opt)
  unique(opt$spcode)
  
  fern25_opt <- fern25_sel[,which(colnames(fern25_sel) %in% opt$spcode==TRUE)]
  colnames(fern25_opt)
  length(fern25_opt[,1])
  
  paFerns <- decostand(fern25_opt[-exx,], method = "pa",1)
  #paFerns <- cbind(fern25[-exx,c(1:2)],paFerns)
  head(paFerns)
  length(paFerns[,1])
  rowSums(paFerns)
  
  opt_mean <- c()
  for (i in seq(paFerns[,1])){
    
    if(rowSums(paFerns[i,])==0) {
      opt_mean[i] <- 0
    } else {
        
    sp <- colnames(paFerns[which(paFerns[i,]!=0)])
    temp_opt <- c()
    for(j in seq(sp))  {
      temp_opt[j] <-  opt[opt$spcode==sp[j],"Optima"]
    }
    opt_mean[i] <- mean(temp_opt)
  }
  }
  
  optSB <- cbind(fern25_sel[-exx,c(1:2)],opt_mean)
  head(optSB)
  length(optSB[,1])
  
  # end of table preparation
  
  # now correlate hand and optSB
  #exclude the hand plots without topographic matchin
  
  optSB <- optSB[optSB$TrNumber %in% trvector,]
  
  toplot2 <-  cbind(toplot,optSB$opt_mean)
  colnames(toplot2)
  
  sbplot1 <- ggplot(toplot2,aes(hand,opt_mean))+
    geom_point(aes(size=3,color=surface))+
    geom_smooth(aes(hand,opt_mean,color=surface),method = "lm",se=FALSE)+
    facet_wrap(~TrNumber, scales="free") +
    theme(strip.text.x = element_text(size=20))+
    theme(strip.text = element_text(size=25))
  sbplot1
  
  sbplot2 <- ggplot(toplot2,aes(hand,opt_mean))+
    geom_point(aes(size=3,color=surface))+
    geom_smooth(aes(hand,opt_mean,color=surface),method = "lm",se=FALSE)+
    theme(strip.text.x = element_text(size=20))+
    theme(strip.text = element_text(size=25))
  sbplot2
 
  
  
