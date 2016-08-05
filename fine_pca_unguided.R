

############# Set pats
ebiopath<-'/ebio'
# ebiopath<-'/home/moisesexpositoalonso/ebio_remote'
ebiopath<-'/Users/moisesexpositoalonso/work_remote/home/moisesexpositoalonso/ebio_remote'

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

acc<-read.csv(paste(ebiopath,"/abt6_projects8/ath_1001G_history/p1001_MOTHER_DATABASE_DOWNSAMPLING/1001G_downsample_genetics_762.csv",sep=""),header=T)

###

dm<-chunkcount
rownames(dm)=dm$Recipient
dm<-dm[,-1]
dm<-apply(dm,2,FUN=function(x){log10(x+1) })

table(dm == Inf)
table(is.na(dm))

pc<-prcomp(dm)
summary(pc)
pc$x
library(devtools)
install_github("ggbiplot", "vqv")
library(ggbiplot)

source("../colors_popgroups.R")

plot(x=pc$x[,1],y=pc$x[,2],col=colors11,pch=16)

g <- ggbiplot(pc, obs.scale = 1, var.scale = 1,groups = acc$kgroup,colour=acc$kgroup, ellipse = TRUE,circle = TRUE)
g <- g + scale_color_manual(colors11)
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')+theme_bw()
