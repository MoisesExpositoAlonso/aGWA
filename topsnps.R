args<-commandArgs(TRUE)

topwanted<-args[1]
#topwanted<-0.0001
#topwanted<-0.1

cat("\n-------------------------------------------------------\n")
cat(paste("\nusage:","Rscript topsnps.R [topwanted] \n"))
cat("\n-------------------------------------------------------\n")


load("results/empresaGWA.RObject")


source('../ancestrygwa_functions.R')
topagwa(empresaGWA,threshold=topwanted)
# topagwa(empresaGWA,threshold=0.0001)
# topagwa(empresaGWA,threshold=100)
# topagwa(empresaGWA,threshold=50)
# topagwa(empresaGWA,threshold=20)
# topagwa(empresaGWA,threshold=10)
# topagwa(empresaGWA,threshold=5)

