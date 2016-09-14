

cat("\n-------------------------------------------------------\n")

cat(paste("\nusage:","Rscript ancestrygwa.R [labelinputfile] [chromosome] [discrete/continuous] "))

cat("\n-------------------------------------------------------\n")

#------------------------------------------------------------------------------------------------------------------------

source("../ancestrygwa_functions.R")

#------------------------------------------------------------------------------------------------------------------------
# Read inputs

args<-commandArgs(TRUE)


labelinputfile<-args[1] 
# labelinputfile<-"label_input_genomes.tsv"

chr<-args[2] 
# chr<-2

#samplefile<-args[3]; #samplefile<-"chromopainterparsedout/chr2-fs.samples.out_parsedGM"
samplefile<-paste0('chromopainterparsedout/',chr,'_parsedGM')

typeofanalysis<-args[3] 
#typeofanalysis<-"discrete"

#------------------------------------------------------------------------------------------------------------------------

# read painted chromosomes
sampleout<-readsampleout(samplefile)

# read phenotypes

labelinput<-readlabelinput(labelinputfile)

typeofvariable(labelinput[,2])

if(length(sampleout[,1]) != length(labelinput[,1])) {
print ("disagreement between painted dataset and labels of genomes!!!")
print(paste("length of sampleout:",length(sampleout[,1]) ) )
print(paste("length of labelinput:",length(labelinput[,1]) ) )

}    #check

# test analyses


print('start aGWA ... ')
timestart<-timestart()

gwares<-apply(sampleout[,-1],MARGIN = 2,FUN = function(x){agwa_test(xmat=x,typeofanalysis)})

cat(paste(' ...end gwa finestructure \n'))
timeend(timestart)


nameobj<-paste(sep="","results/agwa_chr",chr,".RObject")
print(paste("aGWA results stored in: ", nameobj))
save(file=nameobj,gwares)

