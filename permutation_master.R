

cat("\n-------------------------------------------------------\n")

cat(paste("\nusage:","Rscript permuation_master.R [labelinputfile] [chromosome] [replicatename] [discrete/continuous] \n"))

cat("\n-------------------------------------------------------\n")



#------------------------------------------------------------------------------------------------------------------------

source("../ancestrygwa_functions.R")


#------------------------------------------------------------------------------------------------------------------------
## Read iput

args<-commandArgs(TRUE)

labelinputfile<-args[1] 
#abelinputfile<-"label_input_genomes.tsv"

chr<-args[2] 
#chr<-2

replicate<-args[3]

sizesample<-100

samplefile<-paste0('chromopainterparsedout/',chr,'_parsedGM')

typeofanalysis<-args[4] 
#typeofanalysis<-"discrete"

blocksizes<-paste(samplefile,"_BLOCKSIZES",sep="")

cat(paste0("\narguments interpreted: labelinput file: ",labelinputfile,"; chr: ",chr,"; replicate :",replicate,"\n" ))
#------------------------------------------------------------------------------------------------------------------------

# read block length distribution
#if(is.character(blocksizes)){
blocklen<-as.numeric(as.matrix(read.table(blocksizes,header=F)))
#}

# read painted chromosomes
#if(!('sampleout' %in% ls() ) ) {
sampleout<-readsampleout(samplefile)
#}

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


nameobj<-paste(sep="","results/rep",replicate,"_agwa_permuted_chr",chr,".RObject")
print(paste("aGWA results stored in: ", nameobj))
save(file=nameobj,gwares_permuted)



