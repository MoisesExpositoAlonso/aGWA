system('find . -name "agwa_permuted_chr*.RObject" > permutobjects.txt')
aGWAfiles<-as.character.factor(read.table('permutobjects.txt')[,1] )
aGWAfiles

emp<-c()

for (i in aGWAfiles){

     load(aGWAfiles_parse_raw_files[i])
emp<-c(emp,gwares_permuted)

}

gwares_permuted=emp

save(file="agwa_permuted_all.RObject",gwares_permuted)
