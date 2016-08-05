

chr1pos<-t(as.matrix(read.table('finestructure/1135-imp-1.haplotypes-drought_chrpositions',sep=' ',header=F)))
chr2pos<-t(as.matrix(read.table('finestructure/1135-imp-2.haplotypes-drought_chrpositions',sep=' ',header=F)))
chr3pos<-t(as.matrix(read.table('finestructure/1135-imp-3.haplotypes-drought_chrpositions',sep=' ',header=F)))
chr4pos<-t(as.matrix(read.table('finestructure/1135-imp-4.haplotypes-drought_chrpositions',sep=' ',header=F)))
chr5pos<-t(as.matrix(read.table('finestructure/1135-imp-5.haplotypes-drought_chrpositions',sep=' ',header=F)))

# chrpos<-cbind(chr1pos,chr2pos,chr3pos,chr4pos,chr5pos)
chrpos=list(chr1=chr1pos,chr2=chr2pos,chr3=chr3pos,chr4=chr4pos,chr5=chr5pos)

save(chrpos,file = 'chrpos.RObject')

# print(chrpos[1:5,1:5])






