


############# Set pats
 ebiopath<-'/ebio'
#ebiopath<-'/home/moisesexpositoalonso/ebio_remote'
# ebiopath<-'/Users/moisesexpositoalonso/work_remote/home/moisesexpositoalonso/ebio_remote'

# pathresults<-paste( ebiopath,'/abt6_projects8/ath_1001G_history/finestructure/guideddrought' ,sep="")
pathresults<-paste( ebiopath,'/abt6_projects8/ath_1001G_history/finestructure/guidedrelicts' ,sep="")
pathresults<-paste( ebiopath,'/abt6_projects8/ath_1001G_history/finestructure/' ,sep="")



############labels of the painted chromosomes

# donors<-read.table(paste( ebiopath, '/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/population-input-drought.csv' ,sep=""), header=F)
donors<-read.table(paste( ebiopath, '/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/population-input-relicts.csv' ,sep=""), header=F)

donorlabels<-as.character.factor(subset(donors,donors$V2=='D')$V1)
donorlabels_columns<-paste("X",donorlabels,sep='')

labelinput<-read.table(paste( ebiopath, '/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/label-input-relicts.csv' ,sep=""), header=F)
labelinput$V4<-substring(labelinput$V2, 1, 1)
labelinput$V5<-c(1:length(labelinput$V1))

labelinput$V4=NA
for (i in 1:length(labelinput$V2)){
  if(nchar(as.character.factor(labelinput$V2[i])) != 1){
  labelinput$V4[i] = substr(as.character.factor(labelinput$V2[i]), 1, nchar(as.character.factor(labelinput$V2[i]))-1)
  }else{
  labelinput$V4[i]=as.character.factor(labelinput$V2[i])
  }

}

row.names(labelinput)=labelinput$V5
head(labelinput)


#########################################################
chr=2

mu=paste(sep="",'chr',chr,'.mutationprobs.out')
mutprobs=paste(sep="",'chr',chr,'.EMprobs.out')
mutprobs=paste(sep="",'chr',chr,'.chunkcounts.out')
mutprobs=paste(sep="",'chr',chr,'.chunklengths.out')

mutprobs<-read.table(paste(pathresults,paste(sep="",'chr',chr,'.mutationprobs.out'), sep="/") , header=T)
# emprobs<-read.table(paste(pathresults,paste(sep="",'chr',chr,'.EMprobs.out'), sep="/") , header=T)
chunkcount<-read.table(paste(pathresults,paste(sep="",'chr',chr,'.chunkcounts.out'), sep="/") , header=T)
chunklengths<-read.table(paste(pathresults,paste(sep="",'chr',chr,'.chunklengths.out'), sep="/") , header=T)

colnames(mutprobs)[1]<-'id'
colnames(chunkcount)[1]<-'id'
colnames(chunklengths)[1]<-'id'


##################### 762 ###############################

acc<-read.csv(paste( ebiopath,'/abt6_projects8/ath_1001G_history/p1001_MOTHER_DATABASE_DOWNSAMPLING/1001G_downsample_genetics_762.csv' ,sep=""), header=T)
acc

meracc<-merge(listdata[[z]],acc,by="id")

#########################################################

load("../MASTER_DATA.RData")

source('colors_popgroups.R')

listdata=list(mutprobs=mutprobs,chunkcount=chunkcount,chunklengths=chunklengths)
names(listdata)


matrixglm<-matrix(ncol=2+1,nrow=length(donorlabels_columns)  )
colnames(matrixglm)<- c("variable","slope","pvalue")
rownames(matrixglm)=donorlabels_columns
matrixglm
z=2
# for (z in 1:length(listdata)){

merged<-merge(listdata[[z]],dataMerge,by="id")

names(listdata)[z]

# i=donorlabels_columns[1]

for (i in donorlabels_columns){
  lmod<-lm(merged[,"m1d_polqua"] ~ merged[,i])
  pval=summary(lmod)$coefficients[2,4]
  slope=summary(lmod)$coefficients[2,1]
  varname=names(listdata)[z]

  matrixglm[i,]<-c(varname,slope,pval)
  }
# }

matrixglm

plot(merged[,'m1d_polqua'] ~merged[,'X5c'])
plot(merged[,'m1d_polqua'] ~merged[,'X2'])
plot(merged[,'m1d_polqua'] ~merged[,'X7'])
plot(merged[,'m1d_polqua'] ~merged[,'X9'])
plot(merged[,'m1d_polqua'] ~merged[,'X11'])
plot(merged[,'m1d_polqua'] ~merged[,'X5g'])

write.table(sep="\t",quote = F,x=matrixglm,file=paste("matrixglm_",names(listdata)[z] ,".tsv",sep="") )


#### world map

source('/media/moisesexpositoalonso/Data/GoogleDrive/R/2014-11-4_1001g/ggplot_world_map.R')
names(merged)

p<-createmapbase()
p<-createmapbase()+coord_equal(xlim =  c(-26,+90 ), ylim=c(10,70))
# p+geom_point(data=merged,aes(y=latitude.x,x=longitude.x,col=m1d_polqua)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=merged,aes(y=latitude.x,x=longitude.x,col=m1d_polqua)) +scale_colour_gradientn(colours = rainbow(5))
ggsave(filename = "m1d_polqua_distribution.pdf",width=10,height = 6)

pdf("relict_painting_distribution.pdf",width=10,height = 6)
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X5a)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X5b)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X5c)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X5d)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X5e)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X5f)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X5g)) +scale_colour_gradientn(colours = rainbow(5))
dev.off()

p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X1)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X2)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X3)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X4)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X6)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X7)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X8)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X9)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X10)) +scale_colour_gradientn(colours = rainbow(5))
p+geom_point(data=meracc,aes(y=latitude,x=longitude,col=X11)) +scale_colour_gradientn(colours = rainbow(5)) + geom_point(data=meracc[which(meracc$X11==0),],aes(y=latitude,x=longitude),col="red")


meracc[which(meracc$X11==0),]
