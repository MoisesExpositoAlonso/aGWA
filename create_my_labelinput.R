
accmaster<-read.table("../762_accmaster_env_str_drough.tsv",sep="\t",header = T)

wantedcolumns<-c("id","kgroup","m1d_polqua")
label_input_genomes<-accmaster[,wantedcolumns]
write.table(x = label_input_genomes,file = "label_input_genomes.tsv",sep="\t",quote=F,row.names=F,col.names = F)

#names(accmaster)
#wantedcolumns<-c("id","PC3","m1d_polqua")
#label_input_genomes<-accmaster[,wantedcolumns]
#write.table(x = label_input_genomes,file = "label_input_genomes_continuous.tsv",sep="\t",quote=F,row.names=F,col.names = F)

