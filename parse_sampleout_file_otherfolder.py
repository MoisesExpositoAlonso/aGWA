import sys
import subprocess

## define input ".samples.out" file

try:
 samplefile=sys.argv[1]
except IndexError:
 print "no input file"
 print "run as: python parse_sampleout_file_otherfolder.py examplechrX.samples.out" 
 exit()

name=samplefile.split('/')[-1 ]

#samplefile='chr2.samples.out_header'
#samplefile='/ebio/abt6_projects8/ath_1001G_history/finestructure/guideddrought/chr2.samples.out'

s=open(samplefile,'r')


# define the folder where the output file is going to be. If there is not a second argument, then is assumed to be in the current folder.
try:
 newfolder=sys.argv[2]  
except IndexError:
 newfolder="chromopainterparsedout"

call("mkdir %s" %(newfolder))


# the new file is going the be created from the input file, and stored in the newfolder.
gm=[]
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
				# gm.append(" ".join([hap,i[2:50]])) # if you want a subset of all positions just to create a toy output.
		except ValueError:
	 		pass   
for i in gm:
	# print i[0:10]
	gmfile.write(i)

gmfile.close()
