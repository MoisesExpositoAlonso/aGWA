
args<-commandArgs(TRUE)
print(paste("these are the read arguments",args))
# args<-c(1,2,3,4,5)

# chr1pos<-t(as.matrix(read.table('finestructure/1135-imp-1.haplotypes-drought_chrpositions',sep=' ',header=F)))
# chr2pos<-t(as.matrix(read.table('finestructure/1135-imp-2.haplotypes-drought_chrpositions',sep=' ',header=F)))
# chr3pos<-t(as.matrix(read.table('finestructure/1135-imp-3.haplotypes-drought_chrpositions',sep=' ',header=F)))
# chr4pos<-t(as.matrix(read.table('finestructure/1135-imp-4.haplotypes-drought_chrpositions',sep=' ',header=F)))
# chr5pos<-t(as.matrix(read.table('finestructure/1135-imp-5.haplotypes-drought_chrpositions',sep=' ',header=F)))
# chrpos=list(chr1=chr1pos,chr2=chr2pos,chr3=chr3pos,chr4=chr4pos,chr5=chr5pos)

chrpos=list()

for ( i in 1:length(args) ){
chr=args[i]
print(paste("reading chromosome.",chr))
	tmp<-t(as.matrix(read.table(paste0(chr,'_chrpositions'),sep=' ',header=F)))
	# tmp<-t(as.matrix(read.table(paste0(chr,'_chrpositions_head'),sep=' ',header=F))) #PROFILING!
	chrpos[[chr]]<-c(tmp)
}

# chrpos

save(chrpos,file = 'chrpos.RObject')

