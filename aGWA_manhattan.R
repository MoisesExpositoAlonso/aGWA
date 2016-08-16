
cat("\n-------------------------------------------------------\n")
cat(paste("\nusage:","Rscript aGWA_manhattan.R [nameanalysis] \n"))

cat("\n-------------------------------------------------------\n")


load("resaGWA.RObject")
load("empresaGWA.RObject")

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

plotname=paste("aGWA",'all_',nameanalysis,".pdf",sep="_")
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
