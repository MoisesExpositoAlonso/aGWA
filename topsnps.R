args<-commandArgs(TRUE)

topwanted<-args[1]
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
