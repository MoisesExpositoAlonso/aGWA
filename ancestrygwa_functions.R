print(getwd())
cat("\nloading aGWA functions...\n")



# ----------------------------------------------# ----------------------------------------------#

 # ----------------------------------------------#

merge_agwa_chrpos<-function(index,chrpos){
# the raw pvalues matrix to store info
res<-matrix(ncol=3,nrow=0)
colnames(res)=c("pval","chr","pos")

print("merging aGWA pvalues and chromosome-positions...")

for ( i in 1:dim(index)[1] ){
myfile<-index[i,2]
mychr<-as.numeric(index[i,1])
print (index[i,])

    load(myfile)
    if(grepl("empiricalpval", myfile)){
      print("Reading empirical pvalues...")
      tmp<-data.frame(pval=empiricalpval)
      }else{
      print("Reading normal pvalues...")
      tmp<-data.frame(pval=gwares)
    }
    
    tmp$chr<-mychr

print(paste("length of chrpos:",length(chrpos[[mychr ]] )))
print(paste("length of agwa:",length(tmp$chr) ))

    if(length(chrpos[[mychr ]] ) != length(tmp$chr) ){
      print(paste("length of aGWA p values and SNP positions do not coincide for :", myfile))
      if(( length(chrpos[[mychr ]] ) - length(tmp$chr) ) <100) {
     print("the difference is less than 100 SNPs, so forcing merge...")
      tmp$pos<-chrpos[[mychr ]] [1:length(tmp$chr)]

      }else{    next}
    }else{
    tmp$pos<-chrpos[[mychr]]
    res<-rbind(res,tmp)
    }
}

return(res)
}


# ----------------------------------------------#

parseaGWAfiles<-function(aGWAfiles){

aGWAfiles<-sort(aGWAfiles)
aGWAfiles

aGWAfiles_parse<-sub(pattern = paste("agwa_chr",sep=""), replacement = "",x = aGWAfiles )
aGWAfiles_parse<-sub(pattern = "results/", replacement = "",x = aGWAfiles_parse )
aGWAfiles_parse<-sub(pattern = ".RObject", replacement = "",x = aGWAfiles_parse )

aGWAfiles_parse_emp<-aGWAfiles_parse[grep('empirical',aGWAfiles_parse)]
aGWAfiles_parse_emp<-sub(pattern = "_empiricalpval.RObject", replacement = "",x = aGWAfiles_parse_emp )
aGWAfiles_parse_emp_files<-aGWAfiles[grep('empirical',aGWAfiles_parse)]

aGWAfiles_parse_raw<-aGWAfiles_parse[-grep('empirical',aGWAfiles_parse)]
aGWAfiles_parse_raw_files<-aGWAfiles[-grep('empirical',aGWAfiles_parse)]

agwaindex<-cbind(aGWAfiles_parse_raw,aGWAfiles_parse_raw_files)
empresagwaindex<-cbind(aGWAfiles_parse_emp,aGWAfiles_parse_emp_files)

indexes<-list(agwaindex,empresagwaindex)
return(indexes)
}



# ----------------------------------------------

paintedhit<-function(labelinputfile,samplefile,chr,mytop,window=500,nameplot='',scalingfactor=100){
library(reshape2)

  # check
if(mytop$chr != chr){print("the SNP hit table mytop does not coincide with the chr number provided!")
  break
  }

# read painted chromosomes
samplefile<-paste0('chromopainterparsedout/',chr,'_parsedGM')

sampleout<-readsampleout(samplefile)

# read phenotypes

labelinput<-readlabelinput(labelinputfile)

# read pvalues

load(paste0('results/empresaGWA.RObject'))

subp<-empresaGWA[which(as.numeric(empresaGWA$chr) == chr),]
dim(subp)
dim(empresaGWA)
rownames(subp)=1:length(rownames(subp))

# read top hits

cat('
The snp table mytop should look like this:
     SNP  pval chr    pos
2_381223 4e-04   2 381223
')

cat('\n reading and prining plot for:' ,paste(colnames(mytop), mytop,sep=": " ), "\n" )

# find postion in sampleout

focal<-which(subp$pos ==mytop$pos)

sampleoutclean<-sampleout[,-1] # because the first column are the genome names
subsampleout<-sampleoutclean[,c( (focal-window) : (focal+window) )]


# sort the genomes baseed on drought
# sampleout_rank<-subsampleout[order(labelinput[,3]), ]
#sorted based on cluster
sampleout_rank<-subsampleout[order(labelinput[,2]), ]
print("visualization of the painted genome matrix...")
print(sampleout_rank[1:10,1:10])

rownames(sampleout_rank)=c(1:length(rownames(sampleout_rank)))
dim(sampleout_rank)

# paint the matrix


sampleout_rank_painted=t(apply(sampleout_rank,1,FUN=function(x){labelinput[as.numeric(x) , 2]} )  )
sampleout_rank_painted[1:10,1:10]

sampleout_rank_painted.df.melt=melt(sampleout_rank_painted)
head(sampleout_rank_painted.df.melt)
tail(sampleout_rank_painted.df.melt)
colnames(sampleout_rank_painted.df.melt)=c('X1','X2','value')

# sampleout_rank_painted.df.melt$value=as.character.factor(sampleout_rank_painted.df.melt$value)
sampleout_rank_painted.df.melt$X2=sampleout_rank_painted.df.melt$X2-(window+1)
sampleout_rank_painted.df.melt$X1=-sampleout_rank_painted.df.melt$X1
head(sampleout_rank_painted.df.melt)
tail(sampleout_rank_painted.df.melt)

sampleout_rank_painted.df.melt$value<-factor(sampleout_rank_painted.df.melt$value)

# sort snps pvalues

df_sub<-subp[c( (focal-window) : (focal+window) ) , ]
df_sub$cumulative=(1:length(df_sub$pos)) -window+1

## plot

library(ggplot2)
p<-ggplot()+ geom_tile(data = sampleout_rank_painted.df.melt, aes(y=X1, x=X2, fill = value))+ scale_fill_manual(values = mycolors) + theme_bw() +xlab("# SNP position from hit ")+ ylab("genomes sorted by k group     -log10( p-value )") #+ ylim(c(-211, +6*10))
nicebreaks<-round( seq(-dim(sampleout)[1], max(-log10(df_sub$pval)*scalingfactor),by = scalingfactor)/scalingfactor ,digits=0)*scalingfactor

p<-p+geom_point(data = df_sub,aes(x = cumulative,y=-log10(pval)*scalingfactor))  + scale_y_continuous(breaks = nicebreaks,labels =nicebreaks )
p

paintedhitfile=paste("plots/paintedhit/",nameplot,'agwahit_',mytop$chr,"_",mytop$pos,".pdf",sep="")
ggsave(paintedhitfile,width = 7,height = 6)
print(paste("aGWA hit and painted genome stored in: ", paintedhitfile))

}

# ----------------------------------------------


manhattan<-function(data,chrpalette=c('black','black','black','black',"black"),empirical=NULL,type='perchromosome',nameplot='',pdfname=paste0('plots/',nameplot, "aGWA",'manhattan',".pdf")){


print(paste0('printing manhattan plot to ...', pdfname))

if(type =='perchromosome'){
print(paste0('you can also output the manhattan in a single plot with the flag= type="cumulative" '))

pdf(file = pdfname,width=10,height = 2*length(unique(data$chr)) ,useDingbats = F)
par(mfrow=c(length(unique(data$chr)) ,1),mar = c(4,4,1.5,1))

for (chromnum in sort(unique(data$chr) )){

subdata<-subset(data,data$chr==chromnum)
subdata$chr<-as.numeric(subdata$chr)
maxchrompos<-max(data$pos)
#chcol<-c('firebrick3','lightskyblue','palegreen3','orchid4',"tan2")


sizes=-log10(subdata$pval)/max(-log10(subdata$pval))
plot(y=-log10(subdata$pval), x= (subdata$pos / 1e6) ,pch=16,cex=sizes,col=chrpalette[subdata$chr],xlab="Position (Mb)",  ylab= "log10 p-values"  ,xaxt="n",xlim=c(0,maxchrompos/1e6),frame.plot = F)

if(!is.null(empirical) ){
  abline(h=-log10(0.05),col="black",lty=2)
  abline(h=-log10(0.01),col="black",lty=2)
  abline(h=-log10(0.001),col="black",lty=2)
  }else{  abline(h=-log10(0.05/length(data$pval)),col="black",lty=2)}

minc=min(round(subdata$pos/ 1e6))
maxc=max(round(subdata$pos/ 1e6))
seqminmax=seq(minc,maxc,by=1)
axis=axis(side=1, at=seqminmax, labels=seqminmax)
}

dev.off()

}else if(type=="cumulative"){
print(paste0('you can also output the manhattan in a different rows pwer chromosome with the flag= type="perchromosome" '))

pdf(file = pdfname,width=12,height = 5 ,useDingbats = F)

data$sumcum<-1:dim(data)[1]
data$chr<-as.numeric(data$chr)

sizes=-log10(data$pval)/max(-log10(data$pval))
plot(y=-log10(data$pval), x= (data$sumcum ) ,pch=16,cex=sizes,col=chrpalette[data$chr],xlab="Cummulative Position (#SNP)",  ylab= "log10 p-values"  ,xaxt="n",frame.plot = F)

if(!is.null(empirical) ){
  abline(h=-log10(0.05),col="black",lty=2)
  abline(h=-log10(0.01),col="black",lty=2)
  abline(h=-log10(0.001),col="black",lty=2)
  }else{  abline(h=-log10(0.05/length(data$pval)),col="black",lty=2)}

minc=min(round(data$sumcum))
maxc=max(round(data$sumcum))
seqminmax=round(seq(minc,maxc,length.out = 10),digits =0 )
axis=axis(side=1, at=seqminmax, labels=seqminmax)

dev.off()
}
}

# ----------------------------------------------

topagwa<-function(data,threshold, nametable='',tsvname=paste0('tables/',nametable, "toptable_",threshold,".tsv")){

  toptable<-matrix(ncol = dim(data)[2], nrow=0)

for(chromnum in sort(unique(data$chr)) ){

subdata<-subset(data,data$chr==chromnum)
subdata$chr<-as.numeric(subdata$chr)

library(reshape)
df.sorted<-sort_df(subdata,vars='pval')

if(threshold > 1){
  toappend<-head(df.sorted,n=threshold)
    if(any(dim(toappend) ==0)){
      next
    }else{
    toappend$chr=chromnum
    toptable<-rbind(toptable,toappend)
    }
}else if(threshold<1) {
  toappend<-subset(df.sorted,df.sorted$pval<threshold)
    if(any(dim(toappend) ==0)){
      next
    }else{
    toappend$chr=chromnum
    toptable<-rbind(toptable,toappend)
    }
}

} # chromosome loop

if(any(dim(toptable))==0){
  print("NO SNP FOUND PASSING THE THRESHOLD!")
  break
}else{
SNP=paste(toptable$chr,toptable$pos,sep="_")
toptable<-cbind(SNP,toptable)

colnames(toptable)=c("SNP",'pval','chr','pos')

write.table(toptable,file=tsvname,sep='\t',quote = F,row.names = F,col.names = T)

return(toptable)
}
}


# ----------------------------------------------

findroot<-function(file){
  rootfile<-tail( strsplit(file,split = '/')[[1]] ,n = 1)
  return(rootfile)
  }

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

# ----------------------------------------------
timestart<-function(){
ptm <- proc.time()
return(ptm[1])
}
# ----------------------------------------------
timeend<-function(timestart){
  cat(paste('time in sec', c(proc.time()[1] - timestart),
            '\ntime in min', c(proc.time()[1] - timestart)/60,'\n'  ))
}
# ----------------------------------------------
typeofvariable<-function(vec){

print(paste("is numeric",is.numeric(labelinput[,2])))
print(paste("is interger", is.integer(labelinput[,2])))
print(paste("is character",is.character(labelinput[,2])))

return("")
}
# ----------------------------------------------
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
# ----------------------------------------------
readsampleout<-function(samplefile){
cat("\nreading the parsed .sample.out file...\n")

if( paste(sep="",samplefile ,".RObject") %in% paste(sep="/","chromopainterparsedout",list.files(path="chromopainterparsedout")) ) {
cat("\n found the corresponding sample.out.RObject file...\n")

load(paste(sep="",samplefile ,".RObject"))

}else{
cat("\n as matrix ... \n")

sampleout<-read.table(samplefile , header=F,fill=T)
save(sampleout,file=paste(sep="",samplefile ,".RObject"))
}
#cat("this is an example of how it looks a painted chromosome \n")
#print(sampleout[1:10,1:10])
return(sampleout)
cat("\n...done\n")
}
# ----------------------------------------------
agwa_testold<-function(xmat,typeofanalysis){

xmat_translate=labelinput[xmat,2] # important this is translating the painted position with the structure information
paintedSNP=xmat_translate

datanew<-cbind(labelinput,paintedSNP)

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
else if(spearman$estimate <0){  PvalNP= + spearman$p.value} ### change the + for a - to implement directional tests. Be careful because negative p-value will be transformed to NA during log transform and will produce errors. This feature needs quite some debugging

}
return(PvalNP)
# return(Pval)
# return(list(Pval=Pval,Fval=Fval))
}

# ----------------------------------------------

agwa_test<-function(xmat,typeofanalysis,returnrsquare=NULL){

xmat_translate=labelinput[xmat,2] # important this is translating the painted position with the structure information
paintedSNP=xmat_translate

datanew<-cbind(labelinput,paintedSNP)

if(typeofanalysis =="discrete"){

if(length(levels(factor(datanew[,4])) ) ==1){
PvalNP=0
}
else{
# ktest<-kruskal.test(datanew[,3] ~ factor(datanew[,4]) ) # old implementation of non parametric tests
test = lm(datanew[,3] ~ factor(datanew[,4]) )

mypval<-as.matrix(anova(test))[1,5]
myrsq<-summary(test)$adj.r.squared
  
# PvalNP=mypval

}
}
else if(typeofanalysis=="continuous"){

# spearman<-cor.test(x=datanew[,3],y= datanew[,4] ,method="spearman",alternative = "two.sided") # old implementation of non parametric tests

test = lm(datanew[,3] ~ datanew[,4] )

mypval<-as.matrix(anova(test))[1,5]
myrsq<-summary(test)$adj.r.squared
  
# PvalNP=mypval
# if(spearman$estimate >0){ PvalNP=spearman$p.value}
# else if(spearman$estimate <0){  PvalNP= - spearman$p.value} ### change the + for a - to implement directional tests. Be careful because negative p-value will be transformed to NA during log transform and will produce errors. This feature needs quite some debugging

}
# return(PvalNP)
# return(Pval)
# return(list(Pval=Pval,Fval=Fval))


if(is.null(returnrsquare)){
 return(mypval)
} else{ 
  return(list(pval=mypval,r2=myrsq)) 
  }

}

# ----------------------------------------------
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

# ----------------------------------------------
permparallel<-function(vec=c(1:10)){
# ----------------------------------------------
gwares_permuted<-sapply(vec,FUN=function(x){agwa_perm(i=x) })

return(gwares_permuted)
}


# ----------------------------------------------# ----------------------------------------------#
print("Existing functions:")
print(ls())
cat("\n ... done loading functions \n")
