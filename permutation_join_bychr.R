
args<-commandArgs(TRUE)
# args<-c(1,2,3,4,5)

system('find results -name "*agwa_permuted_chr*.RObject" > permutobjects.txt')
aGWAfiles<-as.character.factor(read.table('permutobjects.txt')[,1] )
aGWAfiles



for ( c in 1:length(args) ){
chr=args[c]

filessub<-aGWAfiles[grep(paste0('chr',chr), aGWAfiles)]

emp<-c()

for (i in 1:length(filessub) ){

load(aGWAfiles[i])
emp<-c(emp,gwares_permuted)
}

gwares_permuted=emp

save(file=paste("results/agwa_permuted_chr",chr,".RObject",sep=""),gwares_permuted)

}
