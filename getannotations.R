
args<-commandArgs(TRUE)

toptable<-args[1]
#toptable='toptable_5.tsv'
#toptable='toptable_10.tsv'
#toptable='toptable_20.tsv'
#toptable='toptable_50.tsv'
#toptable='toptable_100.tsv'
toptable='toptable_0.001.tsv'

cat("\n-------------------------------------------------------\n")
cat(paste("\nusage:","Rscript getannotations.R [toptable] \n"))
cat("For Arabidopsis you need these two libraries:
org.At.tair.db
DBI
")
cat("\n-------------------------------------------------------\n")


## get AT annotations

# genelist<-read.table(paste0("tables/",toptable),sep='\t',header=T)
genelist<-read.table(paste0("",toptable),sep='\t',header=T)
genelist=data.frame(chr=genelist$chr,pos=genelist$pos)

source('getannotations_functions.R')
source('../getannotations_functions.R')
newgenelist<-findgenename(snptable = genelist,ebioroot = '~/ebio/')
newgenelist$genenames
cat(unique(newgenelist$genenames),sep="\n")

head(newgenelist)
write.table(sep="\t", x = newgenelist,file=paste0("",toptable,"_newgenelist.tsv"),col.names = T)
write.table(sep="\t", x = na.omit(unique(newgenelist[,3])),file=paste0("",toptable,"_newgenelist_unique.tsv"),col.names = F,row.names = F)

### the human name

library(org.At.tair.db)
library(DBI)
columns(org.At.tair.db)
keytypes(org.At.tair.db)


geneinfolist<-select(org.At.tair.db, keys=na.omit(newgenelist$genenames),columns=c("SYMBOL","GENENAME",'GOALL'))
head(geneinfolist)
geneinfolist<-geneinfolist[c('TAIR', 'SYMBOL' ,'GENENAME','GOALL')]
geneinfolist<-geneinfolist[c('TAIR', 'SYMBOL' )]

cat(unique(geneinfolist$SYMBOL),sep="\n")

geneinfolist_tmp<-as.matrix(tapply(geneinfolist$SYMBOL,geneinfolist$TAIR,FUN = function(x){paste(x,collapse = " | ")}))
newgeneinfolist<-data.frame(TAIR=row.names(geneinfolist_tmp),SYMBOL=geneinfolist_tmp[,1])
newgeneinfolist


write.table(sep="\t", x = newgeneinfolist,file=paste0("",toptable,"_newgeneinfolist.tsv"))

