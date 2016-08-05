import sys

samplefile=sys.argv[1]
name=samplefile.split('/')[-1]
newfolder="/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/unguided"
#samplefile='/ebio/abt6_projects8/ath_1001G_history/finestructure/guideddrought/chr2.samples.out'
# samplefile='chr2.samples.out_header'

s=open(samplefile,'r')

gm=[]
#gmfile=samplefile+'_parsedGM'
gmfile=newfolder+"/"+name+'_parsedGM'
gmfile=open(gmfile,'w')

for i in s:
	if 'HAP' in i:
		split=i.replace('\n','').split(' ')

		hap=split[2]
	else:
		try:
	  		linecode=int(i[0])
	  		if int(i[0]) == 2:
				gm.append(" ".join([hap,i[2:]]))
				# gm.append(" ".join([hap,i[2:50]])) # if you want a subset of all positions
		except ValueError:
	 		pass   

for i in gm:
	# print i[0:10]
	gmfile.write(i)

gmfile.close()


