
cat("\n-------------------------------------------------------\n")
cat(paste("\nusage:","Rscript empirical_distribution_bypermutation.R [labelinputfile] [chromosome] [sampleoutfile] [discrete/continuous] [meanlengthblocks/fileblocklengths] [path_samleout_parsed] [nameanalysis] \n"))
cat(paste("\nrequired fields:","Rscript empirical_distribution_bypermutation.R [compulsory] [compulsory] [sompulsory] [default=discrete] [meanlengthblocks=1000] [optional] [optional] \n"))
cat("\n-------------------------------------------------------\n")

source("ancestrygwa_functions.R")

args<-commandArgs(TRUE)

labelinputfile<-args[1] ; labelinputfile<-"label_input_genomes.tsv"
chr<-args[2] ; #chr<-2
# samplefile<-args[3]; #samplefile<-"chromopainterparsedout/chr2-fs.samples.out_parsedGM"
samplefile<-paste('chromopainterparsedout/chr',chr,'-fs.samples.out_parsedGM',sep='')
typeofanalysis<-args[4] ; typeofanalysis<-"discrete"
nameanalysis<-args[5] ; nameanalysis<-""

#blocksizes<-args[4];blocksizes<-"chromopainterparsedout/chr2-fs.samples.out_parsedGM_BLOCKSIZES"
blocksizes<-paste(samplefile,"_BLOCKSIZES",sep="")


#------------------------------------------------------------------------------------------------------------------------

# read block length distribution
if(is.character(blocksizes)){
blocklen<-as.numeric(as.matrix(read.table(blocksizes,header=F)))
}

# read painted chromosomes
if(!('sampleout' %in% ls() )  )
{
sampleout<-readsampleout(samplefile)
}

# read phenotypes

labelinput<-readlabelinput(labelinputfile)

typeofvariable(labelinput[,2])

if(length(sampleout[,1]) != length(labelinput[,1])) { print ("disagreement between painted dataset and labels of genomes!!!")}    #check

# test analyses

print('start aGWA with permutation ... ')
timestart<-timestart()

# bigass analysis
gwares_permuted<-permparallel(c(1:dim(sampleout)[2] ))
# gwares_permuted<-permparallel(c(1:1000))

cat(paste(' ...end aGWA with permutation. \n'))
timeend(timestart)


nameobj<-paste(sep="",nameanalysis,"agwa_permuted_chr",chr,".RObject")
print(paste("aGWA results stored in: ", nameobj))
save(file=nameobj,gwares_permuted)



