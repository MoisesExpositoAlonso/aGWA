cat(paste("\nusage:","Rscript ancestrygwa.R [labelinputfile] [chromosome] [sampleoutfile] [discrete/continuous] [path_samleout_parsed] [nameanalysis] "))
cat(paste("\nrequired fields:","Rscript ancestrygwa.R [compulsory] [compulsory] [compulsory] [default=discrete] [optional] [optional]"))
cat("\n-------------------------------------------------------\n")

source("ancestrygwa_functions.R")

args<-commandArgs(TRUE)

labelinputfile<-args[1] ; labelinputfile<-"label_input_genomes.tsv"
chr<-args[2] ;# chr<-2
#samplefile<-args[3]; #samplefile<-"chromopainterparsedout/chr2-fs.samples.out_parsedGM"
samplefile<-paste('chromopainterparsedout/chr',chr,'-fs.samples.out_parsedGM',sep='')
typeofanalysis<-args[4] ; typeofanalysis<-"discrete"
nameanalysis<-args[5] ; nameanalysis<-""
#pathresults<-args[6]; pathresults<-"chromopainterparsedout"

#------------------------------------------------------------------------------------------------------------------------

# read painted chromosomes
sampleout<-readsampleout(samplefile)

# read phenotypes

labelinput<-readlabelinput(labelinputfile)

typeofvariable(labelinput[,2])

if(length(sampleout[,1]) != length(labelinput[,1])) { print ("disagreement between painted dataset and labels of genomes!!!")}    #check

# test analyses


print('start aGWA ... ')
timestart<-timestart()

gwares<-apply(sampleout[,-1],MARGIN = 2,FUN = function(x){agwa_test(xmat=x,typeofanalysis)})

cat(paste(' ...end gwa finestructure \n'))
timeend(timestart)


nameobj<-paste(sep="",nameanalysis,"agwa_chr",chr,".RObject")
print(paste("aGWA results stored in: ", nameobj))
save(file=nameobj,gwares)

