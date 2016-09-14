findgenename<-function(snptable=genelist,ebioroot='/home/moisesexpositoalonso/ebio_remote/') {
  
  if(colnames(snptable)[1] != 'chr' & colnames(snptable)[2] != 'pos')  {
    print(' please rename your snptable as two columns: chr pos')
    break
    }
  
  ann<-read.table(paste(ebioroot,'abt6_projects8/ath_1001G_history/TAIR10_GFF3_genes_transposons.gff.txt',sep=''))
  head(ann)
  
  newsnptable<-snptable
  newsnptable$genenames<-NA
for (i in 1:dim(snptable)[1]){
print(i)
  annsub<-subset(ann,ann$V1 == paste("Chr",snptable[i,"chr"],sep='') & ann$V4<snptable[i,"pos"] & ann$V5> snptable[i,"pos"]  )
  annsub_select<-subset(annsub,annsub$V3 %in% c('gene'))
  annsub_select_an<-as.character.factor(annsub_select$V9)
  if(length(annsub_select_an) != 0L){
  genename<-strsplit(annsub_select_an,split = 'Name=')[[1]][2]
  # annsub_collapse<-paste(annsub$V9,collapse = "|")
  # annsub_collapse_1<-gsub(annsub_collapse,pattern = "Name=",replacement = "")
  # annsub_collapse_2<-gsub(annsub_collapse_1,pattern = "ID=",replacement = "")
  # annsub_collapse_3<-gsub(annsub_collapse_2,pattern = "Parent=",replacement = "")
  # strsplit(annsub_collapse_3,split = "|",fix=T)[[1]]
  # grep(pattern = "AT[1-5]G*",x = strsplit( annsub_collapse,split = c(";","=","|") ) )
  # genename<-as.character.factor(annsub$V9[grep("ID=",annsub$V9)] )
  # genename_list<-strsplit(genename,split = ";",fixed = T)[[1]]
  # genename_list_parse<-genename_list[grep("ID=",genename_list)]
  # genename_list_parse_clean<-sub(x = genename_list_parse,pattern = "ID=",replacement = "")
  # newsnptable$genenames[i]<-genename_list_parse
  newsnptable$genenames[i]<-genename
  }else{print('not found')}
  
  }
  return(newsnptable)
}
