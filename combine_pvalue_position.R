## include mapping positions
chrposfile<-list.files(".",pattern="chrpos.RObject")

if(identical(chrposfile, character(0))  ){
cat("problem finding the genome map")
cat("run ")

}else  {
  cat("loading chrpos.RObject file...")
  load(file = chrposfile)
  cat("...done")
}

# head(chrpos)
# tail(chrpos)
# chrpos[["chr1"]]

## aGWA results
nameanalysis<-""

system('find . -name "agwa_chr*.RObject" > agwaobjects.txt')
aGWAfiles<-as.character.factor(read.table('agwaobjects.txt')[,1] )
aGWAfiles

# aGWAfiles<-list.files(path = ".",pattern=paste(nameanalysis,"agwa_chr[[:xdigit:]].RObject",sep=""))
# aGWAfiles

aGWAfiles_parse<-sub(pattern = paste(nameanalysis,"agwa_chr",sep=""), replacement = "",x = aGWAfiles )
aGWAfiles_parse<-sub(pattern = "./", replacement = "",x = aGWAfiles_parse )
aGWAfiles_parse<-sub(pattern = ".RObject", replacement = "",x = aGWAfiles_parse )

aGWAfiles_parse_emp<-aGWAfiles_parse[grep('empirical',aGWAfiles_parse)]
aGWAfiles_parse_emp<-sub(pattern = "_empiricalpval.RObject", replacement = "",x = aGWAfiles_parse_emp )
aGWAfiles_parse_emp_files<-aGWAfiles[grep('empirical',aGWAfiles_parse)]
  
aGWAfiles_parse_raw<-aGWAfiles_parse[-grep('empirical',aGWAfiles_parse)]
aGWAfiles_parse_raw_files<-aGWAfiles[-grep('empirical',aGWAfiles_parse)]

# the raw pvalue
resaGWA<-matrix(ncol=3,nrow=0)
colnames(resaGWA)=c("pval","chr","pos")

for ( i in 1:length(aGWAfiles_parse_raw_files) ){
   load(aGWAfiles_parse_raw_files[i])
# gwares=combinedres # this can be removed when runned again the gwa
    tmp<-data.frame(pval=gwares)
    tmp$chr<-aGWAfiles_parse_raw[i]
    if(length(chrpos[[paste(sep="","chr",aGWAfiles_parse_raw[i])]] ) != length(tmp$chr) ){
      print(paste("length of aGWA p values and SNP positions do not coincide for :", aGWAfiles_parse_raw[i]))
      next
    }
    tmp$pos<-chrpos[[paste(sep="","chr",aGWAfiles_parse_raw[i])]]
    resaGWA<-rbind(resaGWA,tmp)
  }

save(x=resaGWA,file="resaGWA.RObject")

# the empirical pvalue
empresaGWA<-matrix(ncol=3,nrow=0)
colnames(empresaGWA)=c("pval","chr","pos")

for ( i in 1:length(aGWAfiles_parse_emp_files) ){
   load(aGWAfiles_parse_emp_files[i])
    tmp<-data.frame(pval=empiricalpval)
    tmp$chr<-aGWAfiles_parse_emp[i]
    if(length(chrpos[[paste(sep="","chr",aGWAfiles_parse_emp[i])]] ) != length(tmp$chr) ){
      print(paste("length of aGWA p values and SNP positions do not coincide for :", aGWAfiles_parse_emp[i]))
      next
    }
    tmp$pos<-chrpos[[paste(sep="","chr",aGWAfiles_parse_emp[i])]]
    empresaGWA<-rbind(empresaGWA,tmp)
  }

save(x=empresaGWA,file="empresaGWA.RObject")
