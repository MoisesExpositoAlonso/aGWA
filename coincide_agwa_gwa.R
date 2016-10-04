
####### AND THE TOP SNPS FROM PREVIOUSLY

# top<-read.table('../multivargwa/_dLDtopserial/m1d_polqua_LDtop_1_151_gwahits.txt',header = T)
top<-read.table('../multivargwa/_top1000/m1d_polqua_gwahits.txt',header = T)
dim(top)

topagwa<-read.table("tables/toptable_0.001.tsv")
dim(topagwa)

intersect(as.character.factor(top$SNP),as.character.factor(topagwa$V1))

top$SNP<-as.character.factor(top$SNP)
topagwa$V1<-as.character.factor(topagwa$V1)

topagwa$V1round<-unlist(lapply(topagwa$V1,FUN=function(x){ substr(x, 1, nchar(x)-0) }))
top$SNPround<-unlist(lapply(top$SNP,FUN=function(x){ substr(x, 1, nchar(x)-0) }))

intersect(top$SNPround,topagwa$V1round)

genelist=data.frame(chr=c(5,5),pos=c(19029537,19029694))


# 