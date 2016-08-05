

############# Set pats
ebiopath<-'/ebio'
# ebiopath<-'/home/moisesexpositoalonso/ebio_remote'
# ebiopath<-'/Users/moisesexpositoalonso/work_remote/home/moisesexpositoalonso/ebio_remote'

# pathresults<-paste( ebiopath,'/abt6_projects8/ath_1001G_history/finestructure/guideddrought' ,sep="")
pathresults<-paste( ebiopath,'/abt6_projects8/ath_1001G_history/finestructure/weiout/combined/combined_output' ,sep="")


##

namefile="combined_output"

mutprobs<-read.table(paste(pathresults,paste(sep="",namefile,".chunkcounts.out"), sep="/") , header=T)
# emprobs<-read.table(paste(pathresults,paste(sep="",'chr',chr,'.EMprobs.out'), sep="/") , header=T)
mutprobs[1:10,1:10]

chunkcount<-read.table(paste(pathresults,paste(sep="",namefile,'.chunkcounts.out'), sep="/") , header=T)
chunklengths<-read.table(paste(pathresults,paste(sep="",namefile,'.chunklengths.out'), sep="/") , header=T)


### phenotypes

setwd(paste( ebiopath,'/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/' ,sep=""))
load("../MASTER_DATA.RData")

###

dm<-chunkcount
# dm<-chunklengths

dnew<-dm[,-1]
dnew[1:10,1:10]
rownames(dnew)=colnames(dnew)=dm[,1]
dnew[1:10,1:10]
dnew<-as.matrix(dnew)
dim(dnew)
dsub<-dnew[which(rownames(dnew) %in%dataMerge$id ),which(rownames(dnew) %in%dataMerge$id )  ]
dsub[1:10,1:10]
diag(dsub)<-NA

dsubdist<-as.dist(dsub)

dp<-(dist(dataMerge$m1d_polqua,method = 'euclidean'))
dpmat<-as.matrix(dp)
colnames(dpmat)=rownames(dpmat)=dataMerge$id

library(ade4)
mantel.rtest(m1 = dsubdist,m2=dp,nrepet = 500)
plot(pch=16,dpmat~dsub)



hist(dsub)
summary(dsub)

## Correlation of the number of chounks each donor genome gives and the drought resistancve of the reciver
tmpfun<-function(y){
  # ct<-cor.test(method="pearson",x=dataMerge$m1d_polqua, y=log10(y),alternative="greater")
  ct<-cor.test(method="spearman",x=dataMerge$m1d_polqua, y=y,)
  return(ct$p.value)

}
allp<-apply(dsub,2,FUN = function(x){tmpfun(x)})
allp
# tmpfun(dsub[,1])
# ct<-cor.test(method="pearson",x=dataMerge$m1d_polqua, y=log10(dsub[,1]),alternative="greater")

table(allp<0.01)

plot(-log10(allp)~dataMerge$m1d_polqua)
plot(-log10(allp)~dataMerge$kgroup)

dataMerge$allp=allp

library(reshape)
# dataMerge_sort<-sort_df(dataMerge,vars=c("kgroup","allp") )
dataMerge_sort<-sort_df(dataMerge,vars=c("kgroup") )
head(dataMerge_sort)
dataMerge_sort$order<-c(1:dim(dataMerge_sort)[1])

source("../colors_popgroups.R")

pdf("chunkcount_corr_drought_res.pdf",width = 6,height = 5)
boxplot(-log10(allp)~dataMerge$kgroup,col=colors11)
abline(h = -log10(0.05),col="darkgrey",lty=2)
points(-log10(dataMerge_sort$allp)~dataMerge_sort$kgroup,col=colors11[dataMerge_sort$kgroup],pch=16)

plot(-log10(dataMerge_sort$allp)~dataMerge_sort$order,col=colors11[dataMerge_sort$kgroup],pch=16,xlab="genomes",ylab="-log10(p-value)")
abline(h = -log10(0.05),col="darkgrey",lty=2)
dev.off()

cor.test(-log10(allp),dataMerge$m1d_polqua)
plot(-log10(allp)~dataMerge$alt)
plot(-log10(allp)~dataMerge$bio_18)
plot(-log10(allp)~dataMerge$latitude.x)

head(dsub)
hist(log10(dsub))
diag(dsub)<-NA
table(is.na(dsub))

# library(reshape)
# meltdsub<-melt(dsub)
# hist(meltdsub[,3])

min(allp)
which(allp==min(allp))
top<-which(allp<0.005)
top<-names(top)

# toplot<-data.frame(dsub[,"9427"])

pdf("chunkcount_topdrougth_intoothers.pdf",width = 8,height = 4)
par(mfrow=c(1,2))
for (i in 1:length(top)){
  topwhatgroup<-subset(dataMerge,dataMerge$id==top[i])$kgroup
boxplot(dsub[,top[i]] ~ dataMerge$kgroup,col=colors11, main=paste( top[i], "(" , topwhatgroup, ")") ,ylab="chunk count")
boxplot(log10(dsub[,top[i]]) ~ dataMerge$kgroup,col=colors11, main=paste( top[i], "(" , topwhatgroup, ")") ,ylab="log 10 chunk count")
# abline(h = -log10(0.05),col="darkgrey",lty=2)
}
dev.off()


pdf("chunkcount_topdrougth_fromoothers.pdf",width = 8,height = 4)
par(mfrow=c(1,2))
for (i in 1:length(top)){
  topwhatgroup<-subset(dataMerge,dataMerge$id==top[i])$kgroup
boxplot(dsub[top[i],] ~ dataMerge$kgroup,col=colors11, main=paste( top[i], "(" , topwhatgroup, ")") ,ylab="chunk count")
boxplot(log10(dsub[top[i],]) ~ dataMerge$kgroup,col=colors11, main=paste( top[i], "(" , topwhatgroup, ")") ,ylab="log 10 chunk count")
# abline(h = -log10(0.05),col="darkgrey",lty=2)
}
dev.off()



hist(dsub[,"9427"])
which(dsub[,"9427"] > 7 )
topchunks<-names(which(dsub[,"9427"] > 8 ))
plot(dataMerge$m1d_polqua ~dataMerge$T_repro)
points(pch=16,subset(dataMerge,dataMerge$id %in% topchunks)$m1d_polqua~subset(dataMerge,dataMerge$id %in% topchunks)$T_repro)

plot(dataMerge$m1d_polqua ~ log10(dsub[,"9427"]))
plot(x=dataMerge$m1d_polqua ,y=log10(dsub["9427",]))
cor.test(dataMerge$m1d_polqua ,log10(dsub[,"9427"]),alternative="greater")
cor.test(x=dataMerge$m1d_polqua y=log10(dsub["9427",]),alternative="greater")

dsub[,"9427"] == dsub["9427",]



#### Maps of everything

source('/media/moisesexpositoalonso/Data/GoogleDrive/R/2014-11-4_1001g/ggplot_world_map.R')

acc<-read.csv("/media/moisesexpositoalonso/Data/GoogleDrive/INVESTIGACION/PhD/p1001_MOTHER_DATABASE_DOWNSAMPLING/1001G_downsample_genetics_762.csv",header=T)

dsubtomerge<-data.frame(dsub)
dsubtomerge[is.na(dsubtomerge)]<-0
dsubtomerge$id=rownames(dsubtomerge)

merged<-merge(dsubtomerge,acc,by="id")

merged

p<-createmapbase()
p<-createmapbase()+coord_equal(xlim =  c(-26,+90 ), ylim=c(10,70))
# p+geom_point(data=merged,aes(y=latitude.x,x=longitude.x,col=m1d_polqua)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=dataMerge,aes(y=latitude.x,x=longitude.x,col=m1d_polqua)) +scale_colour_gradientn(colours = rainbow(5))


allacc<-colnames(merged)[2:212]

for (i in allacc){

subs=subset(acc,acc$id==gsub(i ,pattern = "X",replacement = "") )
whichacc=subs$kgroup
aci=gsub(i ,pattern = "X",replacement = "")

pdf(paste("all_genomes_painted_locations",paste(whichacc,"-",gsub(i ,pattern = "X",replacement = ""), seo="_"),".pdf"),width=10,height = 6,useDingbats = F)
p+geom_point(data=merged,aes(y=latitude,x=longitude,col=log10(merged[,i]) ))+scale_colour_gradientn(colours = rainbow(5))+ ggtitle(paste(whichacc,"-",gsub(i ,pattern = "X",replacement = ""), sep="_")) + geom_point(data=subs,aes(y=latitude,x=longitude),col="black")
dev.off()

}
