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

############labels of the painted chromoscare in omes


labelinput<-read.table(paste( ebiopath, '/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/label-input-relicts.csv' ,sep=""), header=F)
# labelinput<-read.table(paste( ebiopath, '/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/label-input-drought.csv' ,sep=""), header=F)
labelinput$V4<-substring(labelinput$V2, 1, 1)
labelinput$V5<-c(1:length(labelinput$V1))

labelinput$V4=NA
for (i in 1:length(labelinput$V2)){
  if(nchar(as.character.factor(labelinput$V2[i])) != 1){
  labelinput$V4[i] = substr(as.character.factor(labelinput$V2[i]), 1, nchar(as.character.factor(labelinput$V2[i]))-1)
  }else{
  labelinput$V4[i]=as.character.factor(labelinput$V2[i])
  }

}

row.names(labelinput)=labelinput$V5
head(labelinput)

labelinput$V6<-0;labelinput$V6[which(labelinput$V4 ==11)] <-1 # is north swedish
labelinput$V7<-0;labelinput$V7[which(labelinput$V4 ==9)] <-1 # is mediterranean
labelinput$V8<-0;labelinput$V8[which(labelinput$V4 ==5)] <-1 # is relict
labelinput$V9<-0;labelinput$V9[which(labelinput$V4 ==10)] <-1 # is catalan

######## argument chromosome

args <- commandArgs(trailingOnly = TRUE)
chr=args[1]
# chr=2



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
