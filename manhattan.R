
cat("\n-------------------------------------------------------\n")
cat(paste("\nusage:","Rscript manhattan.R [typeplot] \n"))
cat(" ... by default empirical is plotted \n")
cat("\n-------------------------------------------------------\n")

#------------------------------------------------------------------------------------------------------------------------
args<-commandArgs(TRUE)

typeplot=args[1]
# typeplot<-"perchromosome"
#------------------------------------------------------------------------------------------------------------------------

load("results/resaGWA.RObject")
load("results/empresaGWA.RObject")

#if(dim(resaGWA) == dim(empresaGWA) ){
#  print("careful the dimensions of the raw results and the permutation corrected pvalues do not coincide!")
#}
#

#### PLOT

colorful<-c('firebrick3','lightskyblue','palegreen3','orchid4',"tan2")
onlyblack<-c('black','black','black','black',"black")

source('../ancestrygwa_functions.R')
pdf("plots/relative_pvalue.pdf",width = 5,height = 5,useDingbats = F)
qqGWA(empresaGWA$pval)
dev.off()
manhattan(data=empresaGWA,empirical=T,type=typeplot)
# manhattan(data=empresaGWA,chrpalette=colorful,empirical=T,type='cumulative')
# system("convert -density 300 plot/aGWAmanhattan_perchromosome.pdf plot/aGWAmanhattan_perchromosome.pdf.png")


# trial subsample to see if I can load it in illustrator
source('../ancestrygwa_functions.R')
dim(empresaGWA) # this is very big!
empresaGWA_sample<-empresaGWA[sample(1:dim(empresaGWA)[1],size=100000,prob=-log10(empresaGWA$pval)/max(-log10(empresaGWA$pval)) ) , ] # sampled given the probability of high pvalue, so when you see the plot you don't miss important peaks.
manhattan(data=empresaGWA_sample,empirical=T,type='perchromosome',nameplot="subsampleplot")


pdf("plots/normal_corrected_comparison.R",width = 5,height = 5,useDingbats = F)
plot(-log10(empresaGWA$pval),-log10(resaGWA$pval), ylab="aGWA corrected",xlab="aGWA")
abline(a=0,b=1)
dev.off()
