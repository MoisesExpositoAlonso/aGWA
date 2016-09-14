

############# Set pats
ebiopath<-'~/ebio'

pathresults<-paste( ebiopath,'/abt6_projects8/ath_1001G_history/finestructure/weiout/combined/combined_output' ,sep="")

namefile="combined_output"

############ read chunks

props<-read.table(paste(pathresults,paste(sep="",namefile,".chunkcounts.out"), sep="/") , header=T)
#donors are column names and receivers row names

############ read labels

labelinput<-read.table("label_input_genomes.tsv")

############ prediction of other groups

head(props)[1:10]

prediction_excluding<-function(x,props,labelinput){

whopredicts<-props[,1] [x]
whatgroup<-subset(labelinput,labelinput$V1==whopredicts)[,2]

painting<-props[,x+1]

dm<-cbind(labelinput,painting)

dmsub<-subset(dm,dm[,2] !=whatgroup)

ct<-cor.test(method="spearman",x=dmsub[,3], y=log10(dmsub[,4]) )
toreturn<-as.numeric(ct$p.value * ct$estimate /abs(ct$estimate))
return(toreturn)

}

res<-lapply(1:dim(props[,-1])[2] ,FUN=function(x) prediction_excluding(x,props,labelinput)) # need to remove the column of names

hist(unlist(res))
hist(abs(unlist(res)))

results<-data.frame(labelinput,results=unlist(res))

head(results)

boxplot(abs(results[,4]) ~ results[,2])

source("../colors_popgroups.R")
source("../names_groups.R")

pdf("globalprop.pdf",width = 7, height = 6,useDingbats = F)
boxplot(-log10(abs(results[,4]) )~ results[,2],col=colors11,xaxt='n')
abline(h=-log10(0.05))
axis(1,at = c(1:11),labels = knames,las=2)
dev.off()

# 
# library(ggplot2)
# library(cowplot)
# 
# positiveresults<-results[results[,4] > 0 , ]
# negativeresults<-results[results[,4] < 0 , ]
# 
# 
# ggplot(data=results)+geom_violin(trim=FALSE,aes(y=-log10(abs(results)),x = factor(V2),group=factor(V2),fill=factor(V2)),alpha=0.8)+ theme_cowplot()+ scale_fill_manual(breaks =c(1:11),values=colors11,guide=F ) +xlab("")+ylab("-log10 p-value")
# +geom_boxplot(aes(x = factor(V2),group=V2,y=-log10(abs(results)),fill=V2),color=colors11,width=0.2)
# 
# ggplot(data=positiveresults)+geom_violin(trim=FALSE,aes(y=-log10(abs(results)),x = factor(V2),group=factor(V2),fill=factor(V2)),alpha=0.8)+ theme_cowplot()+ scale_fill_manual(breaks =c(1:11),values=colors11,guide=F ) +xlab("")+ylab("-log10 p-value")
# 
# ggplot(data=negativeresults)+geom_violin(trim=FALSE,aes(y=-log10(abs(results)),x = factor(V2),group=factor(V2),fill=factor(V2)),alpha=0.8)+ theme_cowplot()+ scale_fill_manual(breaks =c(1:11),values=colors11,guide=F ) +xlab("")+ylab("-log10 p-value")


#----------------------------------------------------------------------------------------------------------------------------
## How many chunks individuals from focal group inherit from others??

xchunk=props[1,]
template<-data.frame(donor=c(1:11))

countchunk<-function(x, props, labelinput){
 focal=labelinput[x,] 
 chunkis=props[x,-1]
 newdat<-data.frame(chunkis=as.numeric(chunkis),groups=labelinput[,2])
 newdatsub<-subset(newdat,newdat$groups != focal[,2])
 counting<-tapply(newdatsub$chunkis,newdatsub$groups ,function(x) mean(log10(x+1)) )
 
 countingmelt<-melt(counting)
 colnames(countingmelt) = c("donor","chunks") 
 countingmelt<-merge(template,countingmelt,"donor",all.x=T)
 countingmelt$acc=focal[,1]
countingmelt$receiver=focal[,2]

 countingmelt[is.na(countingmelt$chunks),"chunks"] <-0
 
 return (countingmelt)
 
}

# rescount<-lapply(1:dim(props[,-1])[2])

# rescount<-lapply(1:dim(props[,-1])[2] ,FUN=function(x) countchunk(x,props,labelinput)) # need to remove the column of names

datacount<-matrix(ncol=3,nrow=0)
for (i in 1:dim(props[,-1])[2] ){
  
  vecval<-countchunk(i,props,labelinput)
  datacount=rbind(datacount,vecval)
}

head(datacount)

tail(datacount)

pdf(paste0("prop_donor_for_receiver.pdf"),width = 8, height = 12,useDingbats = F)
par(mfrow=c(4,3))
for(i in 1:11){
# pdf(paste0("prop_donor_for_receiver",i ,".pdf"),width = 7, height = 6,useDingbats = F)

datacountsub<-subset(datacount,datacount$receiver==i)  
boxplot(datacountsub[,2]~ datacountsub[,1],col=colors11,xaxt='n',ylab="log10 # chunks painted",main=knames[i])
# abline(h=-log10(0.05))
# axis(1,at = c(1:11),labels = knames,las=2)
# dev.off()

}

boxplot(-log10(abs(results[,4]) )~ results[,2],col=colors11,xaxt='n')
abline(h=-log10(0.05))
dev.off()
