

cat("\n-------------------------------------------------------\n")

cat(paste("usage:","Rscript combine_pvale_positions.R "))

cat("\n-------------------------------------------------------\n")

#------------------------------------------------------------------------------------------------------------------------

source("../ancestrygwa_functions.R")
# source("ancestrygwa_functions.R")

#------------------------------------------------------------------------------------------------------------------------
## include mapping positions
chrposfile<-paste0("","chrpos.RObject")

if(identical(chrposfile, character(0))  ){
cat("problem finding the genome map")
cat("run ")

}else  {
  cat("loading chrpos.RObject file...")
  load(file = chrposfile)
  cat("...done")
}


## aGWA results

# Parse files
system('find results/ -name "agwa_chr*.RObject" > agwaobjects.txt')
aGWAfiles<-as.character.factor(read.table('agwaobjects.txt')[,1] )

indexes<-parseaGWAfiles(aGWAfiles)

agwaindex<-indexes[[1]]
empresagwaindex<-indexes[[2]]

source("../ancestrygwa_functions.R")

# Merge normal agwa
resaGWA<-merge_agwa_chrpos(index=agwaindex,chrpos)
save(x=resaGWA,file=paste0("results/resaGWA.RObject") )

# Merge empirical pvalue agwa
empresaGWA<-merge_agwa_chrpos(index=empresagwaindex,chrpos)
save(x=empresaGWA,file= paste0("results/empresaGWA.RObject")  )

# 559076
# 314783
# 396501
# 334957
# 475557

print("finished R script for combining pvalues")
