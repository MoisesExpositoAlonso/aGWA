
####### AND THE TOP SNPS FROM PREVIOUSLY

top<-read.table('../multivargwa/_dLDtopserial/m1d_polqua_LDtop_1_151_gwahits.txt',header = T)

topagwa<-read.table("tables/toptable_0.001.tsv")


intersect(as.character.factor(top$SNP),as.character.factor(topagwa$V1))

top$SNP<-as.character.factor(top$SNP)
topagwa$V1<-as.character.factor(topagwa$V1)

topagwa$V1round<-unlist(lapply(topagwa$V1,FUN=function(x){ substr(x, 1, nchar(x)-4) }))
top$SNPround<-unlist(lapply(top$SNP,FUN=function(x){ substr(x, 1, nchar(x)-4) }))

intersect(top$SNPround,topagwa$V1round)
