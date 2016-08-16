
cat("\n-------------------------------------------------------\n")
cat(paste("\nusage:","Rscript empirical_master.R [labelinputfile] [chromosome]  [sizesample=100] [sampleoutfile] [discrete/continuous] [meanlengthblocks/fileblocklengths] [path_samleout_parsed] [nameanalysis/replicate] \n"))

cat("\n-------------------------------------------------------\n")

source("ancestrygwa_functions.R")

args<-commandArgs(TRUE)

labelinputfile<-args[1] ; labelinputfile<-"label_input_genomes.tsv"
chr<-args[2] ; #chr<-2
# samplefile<-args[3]; #samplefile<-"chromopainterparsedout/chr2-fs.samples.out_parsedGM"
sizesample<-args[6]
sizesample<-100

samplefile<-paste('chromopainterparsedout/chr',chr,'-fs.samples.out_parsedGM',sep='')
typeofanalysis<-args[4]
typeofanalysis<-"discrete"

nameanalysis<-args[5]
nameanalysis<-""

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

randomsamplechromosome<-sample(c(1:dim(sampleout)[2]) , size=sizesample)
gwares_permuted<-permparallel( randomsamplechromosome)
# gwares_permuted<-permparallel(c(1:1000))

cat(paste(' ...end aGWA with permutation. \n'))
timeend(timestart)


nameobj<-paste(sep="",nameanalysis,"agwa_permuted_chr",chr,".RObject")
print(paste("aGWA results stored in: ", nameobj))
save(file=nameobj,gwares_permuted)



