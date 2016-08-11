cat("\nloading aGWA functions...\n")


# ----------------------------------------------

qqGWA<-function(pvals,dir=getwd(),name="",save=F){

  observed <- sort(pvals)
  lobs <- -(log10(observed))

  expected <- c(1:length(observed))
  lexp <- -(log10(expected / (length(expected)+1)))

  if(save==T){
    pdf(paste("qqplot",name,".pdf",sep=""), width=6, height=6)
  }

  plot(c(0,7), c(0,7), col="red", lwd=3, type="l", xlab="Expected (-logP)", ylab="Observed (-logP)", xlim=c(0,max(lobs)+1), ylim=c(0,max(lobs)+1), las=1, xaxs="i", yaxs="i", bty="l")
  points(lexp, lobs, pch=23, cex=.4, bg="black")
  abline(lm(lobs~lexp ),col="black",lty="dashed")
  title( paste("lambda",round(coefficients(lm(lobs~lexp ))[2],digits=3) ) )
  if(save==T){
    print(paste("qqplot printed to: ",dir))
    dev.off()
  }}

pvalue_given_distribution<-function(value,distribution){
    # sorted<-sort(distribution)
    pemp<-as.numeric( table(distribution<value)['TRUE'] / length(distribution) )
  if(is.na(pemp)){pemp=0}
  return(pemp)
}


timestart<-function(){
ptm <- proc.time()
return(ptm[1])
}

timeend<-function(timestart){
  cat(paste('time in sec', c(proc.time()[1] - timestart),
            '\ntime in min', c(proc.time()[1] - timestart)/60,'\n'  ))
}

typeofvariable<-function(vec){

print(paste("is numeric",is.numeric(labelinput[,2])))
print(paste("is interger", is.integer(labelinput[,2])))
print(paste("is character",is.character(labelinput[,2])))

return("")
}

readlabelinput<-function(labelinputfile){
  labelinput<-read.table(labelinputfile,header =  F)
}

# readsampleout<-function(pathresults,chr){
  # get the genome matrix files
# listfiles<-list.files(pathresults,pattern="*_parsedGM")
# samplefile<-listfiles[grep(chr,listfiles)]
# cat("found these files: ")
# cat(listfiles)
# cat("using this file: ")
# cat(samplefile)
# # read the sample.out file:
#
# sampleout<-read.table(paste(pathresults,samplefile, sep="/") , header=F,fill=T)
#
# cat("this is an example of how it looks a painted chromosome")
# print(sampleout[1:10,1:10])
# return(sampleout)
# }

readsampleout<-function(samplefile){
cat("\nreading the parsed .sample.out file...\n")
sampleout<-read.table(samplefile , header=F,fill=T)

cat("this is an example of how it looks a painted chromosome \n")
print(sampleout[1:10,1:10])
return(sampleout)
cat("\n...done\n")
}

agwa_test<-function(xmat,typeofanalysis){

xmat_translate=labelinput[xmat,2] # important this is translating the painted position with the structure information
paintedSNP=xmat_translate

datanew<-cbind(labelinput,paintedSNP)

# colanmes(datanew)<-c("id","structure","phenotype" ,"paintedSNP") ### Just as explanation, these are the columns

# ANALYSIS

if(typeofanalysis =="discrete"){

if(length(levels(factor(datanew[,4])) ) ==1){
PvalNP=0
}
else{
ktest<-kruskal.test(datanew[,3] ~ factor(datanew[,4]) )

PvalNP=ktest$p.value
}
}
else if(typeofanalysis=="continuous"){

spearman<-cor.test(x=datanew[,3],y= datanew[,4] ,method="spearman",alternative = "two.sided")

if(spearman$estimate >0){ PvalNP=spearman$p.value}
else if(spearman$estimate <0){  PvalNP= - spearman$p.value}

}

return(PvalNP)
# return(Pval)
# return(list(Pval=Pval,Fval=Fval))

}


#agwa_permut<-function(){
#
#gwares_permuted<-c()
#
#
##for(i in 1:dim(sampleout)[2] ){
#for(i in 1:10 ){ ################ change thisss!!!!!!!!!!!!!!!
#
#if (length(blocklen ) >1)
#{
#  vectorpermutation<-as.numeric(sample(size=dim(sampleout)[1],x =round(as.numeric(names(table(blocklen)))/2,digits = 0),prob = table(blocklen) ,replace = T) )
#}
#else if(is.numeric(blocklen ) & length(blocklen) ==1)
#{
#vectorpermutation<-sample(size=dim(sampleout)[1],x =c(-1,1),prob = c(0.5,0.5),replace = T ) *round(blocklenfa/2,digits=0)
#}
#  else{cat("something wrong with the block length distribution or value!")}
#
##print (vectorpermutation)
#
#positionsinmatrix<- abs(i - vectorpermutation) +1 #just in case sometimes the randomization falls outside the chromosome because is in one extreme. also +1 because the first column are the names of chromosomes, the iteration should start in 2
#
##print(positionsinmatrix)
#
##xmat<-sampleout[1:dim(sampleout)[1],positionsinmatrix]
#xmat<-sampleout[cbind(1:dim(sampleout)[1],positionsinmatrix)]
#
#print(xmat)
#gwares_permuted[i]<-agwa_test(xmat=xmat,typeofanalysis)
#
#print(gwares_permuted)
#
#}
#
##print (head(gwares_permuted))
##print (tail(gwares_permuted))
#
#}

agwa_perm<-function(i){

#if (length(blocklen ) >1)
#{
  vectorpermutation<-as.numeric(sample(size=dim(sampleout)[1],x =round(as.numeric(names(table(blocklen)))/2,digits = 0),prob = table(blocklen) ,replace = T) )
#}
#else if(is.numeric(blocklen ) & length(blocklen) ==1)
#{
#vectorpermutation<-sample(size=dim(sampleout)[1],x =c(-1,1),prob = c(0.5,0.5),replace = T ) *round(blocklenfa/2,digits=0)
#}
#  else{cat("something wrong with the block length distribution or value!")}

#print (vectorpermutation)

positionsinmatrix<- abs(i - vectorpermutation) +2 #just in case sometimes the randomization falls outside the chromosome because is in one extreme. also +1 because the first column are the names of chromosomes, the iteration should start in 2

#print(positionsinmatrix)
#print(any(positionsinmatrix ==1))
#xmat<-sampleout[1:dim(sampleout)[1],positionsinmatrix]

xmat<-sampleout[cbind(1:dim(sampleout)[1],positionsinmatrix)]

#print(xmat)
snpres<-agwa_test(xmat=xmat,typeofanalysis)

return(snpres)
}


permparallel<-function(vec=c(1:10)){

gwares_permuted<-sapply(vec,FUN=function(x){agwa_perm(i=x) })

return(gwares_permuted)
}


# ----------------------------------------------
cat("\n ... done loading functions \n")
