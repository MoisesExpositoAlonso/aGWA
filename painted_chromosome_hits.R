# painted chromosome hits
finetop<-read.table('finegwa_TOP_.csv',sep='\t',header=T)
finetop
library(reshape)
finetop_sort<-sort_df(finetop,vars='fd')
head(finetop_sort,n=50)

i=1
for (i in 1:length(finetop_sort$fd)){
focalsnp<-finetop_sort[i,]

###########
chr=focalsnp$chr


############# Set pats
ebiopath<-'/ebio'
# ebiopath<-'/home/moisesexpositoalonso/ebio_remote'
# ebiopath<-'/Users/moisesexpositoalonso/work_remote/home/moisesexpositoalonso/ebio_remote'

pathresults<-paste( ebiopath,'/abt6_projects8/ath_1001G_history/finestructure/guideddrought' ,sep="")


############labels of the painted chromosomes

labelinput<-read.table(paste( ebiopath, '/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/label-input-drought.csv' ,sep=""), header=F)
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

############ colors per group

source(paste( ebiopath,'/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/ancestrytest/colors_popgroups.R' ,sep="") )

colors11

############ load drought data

load('../MASTER_DATA.RData')

dataMerge
rankdrought<-sort_df(dataMerge,vars='m1d_polqua')$id

sort_df(dataMerge,vars='m1d_polqua')[,c("id","m1d_polqua")]

#### load GENOME MATRIX
samplefile=paste(sep="",'chr',chr,'.samples.out_parsedGM')

# samplefile='chr2.samples.out_parsedGM'
# sampleout<-read.table(paste(pathresults,samplefile, sep="/") , header=F,fill=T)
# print(sampleout[1:10,1:10])

#### subsetted in terminal
#### subset genome matrix
load('chrpos.RObject')
positions=chrpos[[chr]]
whichindex=which(positions == focalsnp$pos)
window=200
# diffposlow=positions-(focalsnp$pos - window)
# diffpowupp=positions-(focalsnp$pos + window)
# whichindex_low<-positions[which(diffposlow<0 &&   == min(abs(positions-focalsnp$pos - window)) )) ]
# whichindex_upp<-positions[which(abs(positions-focalsnp$pos + window) == min(abs(positions-focalsnp$pos + window) ))]
# whichindex_low<focalsnp$pos
# whichindex_low=which(positions == focalsnp$pos - window)
# whichindex_up=which(positions == focalsnp$pos + window)

samplefilepath<-paste(pathresults,samplefile, sep="/")
samplefilepathout<-paste(samplefilepath,"_subsetforplot_",focalsnp$pos,sep="")
comand<-paste('cat ', samplefilepath,' | cut -d " " -f ' ,whichindex-window,"-", whichindex+window , ">" ,samplefilepathout, sep='')
cat(comand)
system(comand)
# cat /ebio/abt6_projects8/ath_1001G_history/finestructure/guideddrought/chr2.samples.out_parsedGM | cut -d " " -f 68159-68359

sampleout<-read.table(samplefilepathout , header=F,fill=T)
dim(sampleout)
rownames(sampleout)=dataMerge$id
print(sampleout[1:10,1:10])

# ranked based in drought

sampleout_rank<-sampleout[order(rankdrought), ] 
print(sampleout_rank[1:10,1:10])
rownames(sampleout_rank)=c(1:length(rownames(sampleout_rank)))
dim(sampleout_rank)

#### plot baseed in color
library('ggplot2')

# sampleout_rank[1,] 
labelinput[as.numeric(sampleout_rank[1,]) , "V4"]
labelinput[as.numeric(sampleout_rank[2,]) , "V4"]

labelinput$V4

sampleout_rank_painted=t(apply(sampleout_rank,1,FUN=function(x){labelinput[as.numeric(x) , "V4"]} )  )
sampleout_rank_painted[1:10,1:10]

sampleout_rank_painted.df.melt=melt(sampleout_rank_painted)
head(sampleout_rank_painted.df.melt)
tail(sampleout_rank_painted.df.melt)
sampleout_rank_painted.df.melt$value=as.character.factor(sampleout_rank_painted.df.melt$value)
sampleout_rank_painted.df.melt$X2=sampleout_rank_painted.df.melt$X2-(window+1)
sampleout_rank_painted.df.melt$X1=-sampleout_rank_painted.df.melt$X1

newcolors<-colors11[sort(unique(as.numeric(sampleout_rank_painted.df.melt$value) ))]


# ggplot(data = sampleout_rank_painted.df.melt, aes(y=X1, x=X2, fill = value))+ geom_tile()+ scale_fill_manual(values = newcolors) + theme_bw()+xlab("# SNP powition from hit ")+ ylab("top - bottom drought resistance accessions ") + ylim(c(-211, +6))
# ggsave(paste('finehit_paintedregion_',focalsnp$chr,"_",focalsnp$pos,".pdf",sep=""),width = 7,height = 6)


# load all finegwas and get the equivalent on top

load('finegwa_chr1.RObject')
f1=res
load('finegwa_chr2.RObject')
f2=res
load('finegwa_chr3.RObject')
f3=res
load('finegwa_chr4.RObject')
f4=res
load('finegwa_chr5.RObject')
f5=res

fd<-list(f1,f2,f3,f4,f5)

df<-data.frame(fd=fd[[focalsnp$chr]],pos=chrpos[[focalsnp$chr]])
df_sub<-df[c((whichindex-window):(whichindex+window)),]
dim(df_sub)
df_sub$cumulative=(1:length(df_sub$pos)) -window+1
  
# ggplot(data = df_sub,aes(x = cumulative,y=-log10(fd)))+ xlab("relative position") + ylab(" - log10 (p value) ")+ geom_point()+theme_bw()
# ggsave(paste('finehit_pvalue_',focalsnp$chr,"_",focalsnp$pos,".pdf",sep=""),width = 7,height = 6)

p<-ggplot()+ geom_tile(data = sampleout_rank_painted.df.melt, aes(y=X1, x=X2, fill = value))+ scale_fill_manual(values = newcolors) + theme_bw()+xlab("# SNP powition from hit ")+ ylab("top - bottom drought resistance accessions       -log10( p-value )") + ylim(c(-211, +6*10))
p+geom_point(data = df_sub,aes(x = cumulative,y=-log10(fd)*10))  + scale_y_continuous(breaks = seq(-210,+60,by=10),labels =c(seq(-210,0,by=10),c(1:6)) )

ggsave(paste('finehit_',focalsnp$chr,"_",focalsnp$pos,".pdf",sep=""),width = 7,height = 6)

}

