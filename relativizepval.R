source('ancestrygwa_functions.R')

args<-commandArgs(TRUE)
chr<-args[1]; chr=2
#agwafile<-args[1] ; agwafile<-'agwa_chr2.RObject'
agwafile<-paste('agwa_chr',chr,'.RObject',sep="")
#permfile<-args[2] ; permfile<-'agwa_permuted_chr2.RObject'
# permfile<-paste( 'agwa_permuted_chr',chr,'.RObject')
permfile<-paste( 'agwa_permuted_all.RObject')


load(agwafile)
agwa=gwares
pdf(paste(agwafile,".hist.pdf",sep=''))
hist(as.numeric(as.matrix(agwa)))
qqGWA(agwa)
dev.off()


load(permfile)
permgwa=gwares_permuted
pdf(paste(permfile,".hist.pdf",sep=''))
hist(as.numeric(as.matrix(permgwa)))
qqGWA(permgwa)
dev.off()

print('start empirical pval distribution correction ... ')
timestart<-timestart()

# empiricalpval<-sapply(agwa[1:10],FUN=function(x){pvalue_given_distribution(x,permgwa) }) # profiling
empiricalpval<-sapply(agwa,FUN=function(x){pvalue_given_distribution(x,permgwa) })
empiricalpval[empiricalpval==0]<- min(empiricalpval[empiricalpval!=0]) # because there can be p values 0 (can mean two things, that the all individuals have an allele state coming from the same population, or that the p value is so low that cannot be achieved by permutation.). In any case I arbitrarily asign those the minimum next pvalue to avoid infinite numbers when transforming

cat(paste(' ...end empirical pval distribution correction. \n'))
timeend(timestart)

empiricalfile<-paste(agwafile,"_empiricalpval",sep="")
save(empiricalpval,file=paste(empiricalfile,".RObject",sep='') )

pdf(file=paste(empiricalfile,".hist.pdf",sep=''))
hist(x = as.numeric(as.matrix(empiricalpval)),col = 'black',breaks=50)
plot(-log10(as.numeric(as.matrix(empiricalpval)),pch=16))
qqGWA(empiricalpval)
dev.off()

