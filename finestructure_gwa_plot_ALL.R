# cent1<-c(1,	13700000,	15900000,	2200001,	30432563)
# cent2<-c(2,	2450000,	5500000,	3050001	,19705359)
# cent3<-c(3,	11300000,	14300000,	3000001	,23470805)
# cent4<-c(4,	1800000,	5150000,	3350001	,18585042)
# cent5<-c(5,	11000000,	13350000,	2350001	,26992728)
# 
# cent<-rbind(cent1,cent2,cent3,cent4,cent5)
getwd()

# args <- commandArgs(trailingOnly = TRUE)
# gwatype=args[1]

# gwatype<-'res'
# gwatype<-'resNSwed'
# gwatype<-'resRel'
# gwatype<-'resMed'
# gwatype<-'resCat'

for (gwatype in c("res","resNSwed","resRel","resMed","resCat")){

print(gwatype)

# load chrpos positions

load('chrpos.RObject')


# load all finegwas
load('finegwa_chr1.RObject')
res=combinedres[,gwatype]
f1=res
load('finegwa_chr2.RObject')
res=combinedres[,gwatype]
f2=res
load('finegwa_chr3.RObject')
res=combinedres[,gwatype]
f3=res
load('finegwa_chr4.RObject')
res=combinedres[,gwatype]
f4=res
load('finegwa_chr5.RObject')
res=combinedres[,gwatype]
f5=res

fd<-list(f1,f2,f3,f4,f5)
head(f1)
head(f2)     
head(f3)
head(f4)
head(f5)

print ("check combinedres")
# head(combinedres)

print ("check fd")
# head(fd)

#length(fd[[chromnum]]) == length(chrpos[[chromnum]])

# setupt plot

Chromosome1<-max(chrpos[[1]])
chromnum=1
col='black'

#chcol<-c('firebrick3','lightskyblue','palegreen3','orchid4',"tan2")
chcol<-c('black','black','black','black',"black")

plotname=paste("fineplot",'all_',gwatype,".pdf",sep="_")
print(plotname)

pdf(file = plotname,width=10,height = 10,useDingbats = F)

par(mfrow=c(5,1),mar = c(4,4,1.5,1))

chromnum=1

for (chromnum in c(1,2,3,4,5)){
  
if(  length(fd[[chromnum]]) != length(chrpos[[chromnum]]) ){ print("SOMETHING WENT WRONG!"); break}

plot(y=-log10(fd[[chromnum]]+1e-7), x= (chrpos[[chromnum]] / 1e6) ,pch=16,cex=-log10(fd[[chromnum]]+1e-7)/max(-log10(fd[[chromnum]]+1e-7)),col=chcol[chromnum],xlab="Position (Mb)",  ylab= "log10 p-values"  ,xaxt="n",xlim=c(0,Chromosome1/1e6),frame.plot = F)
abline(h=-log10(0.05/length(fd[[chromnum]])),col="grey")
minc=min(round(chrpos[[chromnum]]/ 1e6))
maxc=max(round(chrpos[[chromnum]] / 1e6))
seqminmax=seq(minc,maxc,by=1)
axis=axis(side=1, at=seqminmax, labels=seqminmax)

}

dev.off()


topwanted=10

for(chromnum in c(1,2,3,5)){
df<-data.frame(fd=fd[[chromnum]],pos=chrpos[[chromnum]])
library(reshape)
df.sorted<-sort_df(df,vars='fd')
if (chromnum==1){
top100<-head(df.sorted,n=topwanted)
top100$chr=chromnum
}else{
  toappend<-head(df.sorted,n=topwanted)
  toappend$chr=chromnum
  top100<-rbind(top100,toappend)
}
}

SNP=paste(top100$chr,top100$pos,sep="_")
top100<-cbind(SNP,top100)
print (top100)
tablename=paste('finegwa_TOP_',gwatype,'.csv',sep="")
print(tablename)
write.table(top100,file=tablename,sep='\t',quote = F,row.names = F,col.names = T)



}