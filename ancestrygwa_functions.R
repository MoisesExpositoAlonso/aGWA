cat("loading aGWA functions...\n")
# load(file = "chrpos.RObject")

aGWA_test<-function(xmat,var="V4"){

xmat_translate=labelinput[xmat,var]

# subGM<-cbind(sampleout[,1],xmat)
subGM<-cbind(sampleout[,1],xmat_translate)

# subGM<-cbind(labelinput$V1,xmat)  # DO NOT USE, JUST PROFILING BEAUSE I USED A REDUCED SAMPLEOUT EXAMPLE
# subGM<-cbind(labelinput$V1,xmat_translate) # DO NOT USE, JUST PROFILING BEAUSE I USED A REDUCED SAMPLEOUT EXAMPLE

head(subGM);tail(subGM)
colnames(subGM)<-c('V1','V2')

# Generate pheno file
load('../MASTER_DATA.RData')
phenoname="m1d_polqua"
mypheno<-dataMerge[,c('id',phenoname)]


# merge subset for anova analysis

datanew<-merge(mypheno,by.x = 'id',subGM,by.y='V1')

# ANALYSIS

#anov<-lm(datanew[,phenoname] ~ factor(datanew[,3]) )
# print(datanew[,3])
#anovres<-anova(anov)

if(length(levels(factor(datanew[,3])) ) ==1){
PvalNP=0
}
else{
ktest<-kruskal.test(datanew[,phenoname] ~ factor(datanew[,3]) )

PvalNP=ktest$p.value
# print(anovres)
#Pval=anovres$`Pr(>F)` [1]
#Fval=anovres$`F value` [1]

# cat(c('Fval','\t','Pval','\n',Fval,'\t',Pval))
}

return(PvalNP)
# return(Pval)
# return(list(Pval=Pval,Fval=Fval))

}
