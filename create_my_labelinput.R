
accmaster<-read.table("../752_accmaster_env_str_drough.tsv",sep="\t",header = T)

wantedcolumns<-c("id","kgroup","m1d_polqua")

label_input_genomes<-accmaster[,wantedcolumns]

write.table(x = label_input_genomes,file = "label_input_genomes.tsv",sep="\t",quote=F,row.names=F,col.names = F)


