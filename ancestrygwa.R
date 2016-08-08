source("ancestrygwa_functions.R")

args<-commandArgs(TRUE)

labelinputfile<-args[1]
chr<-args[2]

pathresults<-"chromopainterparsedout"

# get the genome matrix files
listfiles<-list.files(pathresults,pattern="_parsedGM")
myfile<-gstfiles[grep(chr,listfiles)]
cat("working over these files:")
cat(listfiles)

# read the sample.out file:

cat("this is an example of how it looks a painted chromosome")
print(sampleout[1:10,1:10])

sampleout<-read.table(paste(pathresults,samplefile, sep="/") , header=F,fill=T)



#### load GENOME MATRIX
# samplefile=paste(sep="",'chr',chr,'.samples.out_parsedGM')
samplefile=paste(sep="",'chr',chr,'-fs.samples.out_parsedGM')
#
# # samplefile='chr2.samples.out_parsedGM'
sampleout<-read.table(paste(pathresults,samplefile, sep="/") , header=F,fill=T)
print(sampleout[1:10,1:10])




# test analyses

print('start gwa finestructure ... ')
res<-apply(sampleout[,-1],MARGIN = 2,FUN = function(x){fineANOVA(xmat=x)})
resNSwed<-apply(sampleout[,-1],MARGIN = 2,FUN = function(x){fineANOVA(xmat=x,var='V6')})
resMed<-apply(sampleout[,-1],MARGIN = 2,FUN = function(x){fineANOVA(xmat=x,var='V7')})
resRel<-apply(sampleout[,-1],MARGIN = 2,FUN = function(x){fineANOVA(xmat=x,var='V8')})
resCat<-apply(sampleout[,-1],MARGIN = 2,FUN = function(x){fineANOVA(xmat=x,var='V9')})

print('end gwa finestructure')

print (head(res))
print (tail(res))

combinedres<-cbind(res,resNSwed,resMed,resRel,resCat)

nameobj<-paste(sep="","finegwa_chr",chr,".RObject")
# save(file=nameobj,res)
print(nameobj)
save(file=nameobj,combinedres)

#
#head(combinedres,n=100)
