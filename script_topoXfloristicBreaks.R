library(mvpart)
head(fern25)
topography <- read.csv("topography.csv", header = TRUE, sep = ",", dec = ".", na.strings = NA)
head(topography)

#transects to exclude in ferns and topography
as.integer(names(which(tapply(fern25$subunit,fern25$TrNumber,length)!=20)))
as.integer(names(which(tapply(topography$topo,topography$tr,length)<11)))
toexclude <- c("776","780","782","794","798","805","806","807","808")


fern25_sel <- fern25[which(! fern25$TrNumber %in% toexclude),]
trsel <- topography[which(! topography$tr %in% toexclude),]


# calculando as distâncias

mfrow <- c(9,9)
pdf(file="speciesBreaks_sub.pdf", height = 7*mfrow[1], width = 8*mfrow[2])
par(mfrow=mfrow, mar=c(3,3,2,0.5), oma=c(2,1,.5,0.5), mgp=c(1.7,0.6,0))

for(k in unique(fern25_sel$TrNumber)){
  
 
k=752
fern25_k <- subset(fern25_sel, fern25_sel$TrNumber==k)
distteste <- vegdist(decostand(fern25_k[-c(1:2)], method = "pa",1), method = "bray"); distteste
distteste <- as.matrix(distteste)

if(F){
dists_extract <- c()
totext <- c()
for(i in seq(1:19)){
  dists_extract[i] <- distteste[i+1,i]
  totext[i] <- paste(i,i+1)
}; dists_extract

}

topo_k <- subset(trsel, trsel$tr==k)
#plot(topo_k$sub,topo_k$topo,type="l",col="red",main=k,ylab = "",yaxt='n',xlab="meters",cex.main=2,cex.lab=2)
#par(new=T)
#plot(dists_extract,type = "h",xaxt='n',ylab = "distance floristic between subunits",xlab="",cex.lab=2 )
#text(dists_extract, y=NULL,totext,cex=2)


#dev.off()

if(F){
plot(topo_k$topography,dists_extract)
cut(topo_k$topography,3)
cut(dists_extract,3)
}

################## trying the MRT as chronological cluster
dist <- c(1:20 * 25)
topoaprox <- approx(topo_k$sub,topo_k$topo,dist)$y
trname <- rep(k,20)
topodata <- cbind(trname,dist,topoaprox)
colnames(topodata)[1:2] <- c("tr","sub")

topodata <- as.data.frame(topodata)
#plot(dist, topoaprox, type="l")

mrttree <- mvpart(distteste~sub,data=topodata,xv="1se",xval = nrow(distteste), xvmult = 100, which=4,margin=0.08,size=4,plot.add = FALSE)
summary(mrttree)
#plot(mrttree)
wheremeter <- data.frame(seq(25,500,25),mrttree$where)
colnames(wheremeter) <- c("sub","class")

plot(topodata$sub, topodata$topoaprox, type="l",ylim=c(0,max(topodata$topoaprox)),main=k,xlab="meters",ylab="topo",cex.lab=2,cex.main=2)
par(new=T)
plot(wheremeter$sub,rep(1,20),col=as.numeric(wheremeter$class),type="p",pch="-",ylab="",cex=5,
     xlab="",ylim=c(0,max(topodata$topoaprox)))

############################################################################
# fim do looping

}

dev.off()



