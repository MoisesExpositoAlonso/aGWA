


cat("\n-------------------------------------------------------\n")

cat(paste("\nusage:","Rscript paintedhits.R [labelinputfile] [discrete/continuous] [window] [toptable]"))

cat("\n-------------------------------------------------------\n")


#------------------------------------------------------------------------------------------------------------------------

source("../ancestrygwa_functions.R")

#------------------------------------------------------------------------------------------------------------------------
# Read inputs

args<-commandArgs(TRUE)


labelinputfile<-args[1] 
#labelinputfile<-"label_input_genomes.tsv"

typeofanalysis<-args[2] 
#typeofanalysis<-"discrete"


window=args[3] 


toptable=args[4]


# -------------------------------------------------------------------------------------------------------


mycolors<-c('#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','gold','#a6cee3')


if( !any(ls()=="mycolors")) {
  if(typeofanalysis=="discrete"){
  mycolors<-rainbow(length(unique(labelinput[,2])))
  print("vector of mycolors not found so generated from rainbow!")
  }else if(typeofanalysis=='continuous'){
  mypallete<-colorRampPalette(c("cyan","navy","black","red","orange"))
  mycolors<-mypallete(dim(labelinput)[1]) [rank(labelinput[,3])]
  }
  else{
    break
    print("type of analysis not provided. Choose either discrete or continuous")
  }
}else{
  print('color pallete found ... did not check the mycolors though')
}
  

# -------------------------------------------------------------------------------------------------------

top<-read.table(file = paste0('',toptable )  ,sep='\t',header=T)
head(top)

# round positions so do I don't repeat plots of very close SNPs (Mbase scale)
top$posround<-round(top$pos/1000)*1000

# remove duplicates
topnoduplicates<-top[!duplicated(top[,"posround"]),]

if(dim(topnoduplicates)[1] > 50 ){
  topnoduplicates<-topnoduplicates[sample(1:dim(topnoduplicates)[1], size=50 ),]
}

# -------------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------------

source("../ancestrygwa_functions.R")

system("mkdir plots/paintedhit/")

for(i in 1:dim(topnoduplicates)[1]){

selected<-topnoduplicates[i,]
paintedhit(labelinputfile,samplefile,chr = selected$chr ,mytop = selected)

}

# which(topnoduplicates$SNP =='5_7575534')
# i=202
