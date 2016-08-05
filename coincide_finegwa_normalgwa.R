# how many are shared


############# Set pats
ebiopath<-'/ebio'
# ebiopath<-'/home/moisesexpositoalonso/ebio_remote'
# ebiopath<-'/Users/moisesexpositoalonso/work_remote/home/moisesexpositoalonso/ebio_remote'

# pathresults<-paste( ebiopath,'/abt6_projects8/ath_1001G_history/finestructure/guideddrought' ,sep="")
pathresults<-paste( ebiopath,'/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/' ,sep="")



fine<-read.table(paste(pathresults,"fineancestry/finegwa_TOP_res.csv",sep=""),header=T)
topgwa<-read.table(paste(pathresults,"polquagwa/m1d_polqua.csv.ps_topgwa_100.tsv",sep=""),header=T)

fine
topgwa

length(setdiff(fine$SNP,topgwa$SNP)) / length(fine$SNP)
