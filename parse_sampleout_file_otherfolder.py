import sys
from subprocess import *

## define input ".samples.out" file

try:
 samplefile=sys.argv[1]
except IndexError:
 print "no input file"
 print "run as: python parse_sampleout_file_otherfolder.py examplechrX.samples.out X" 
 exit()

name=samplefile.split('/')[-1 ]

print samplefile
s=open(samplefile,'r')

## get the chromosome name
try:
 chrom=sys.argv[2]
except IndexError:
 print "no input chromosome name"


#define the folder where the output file is going to be. If there is not a second argument, then is assumed to be in the current folder.
# try:
 # newfolder=sys.argv[2]  
# except IndexError:  # I have fixed the folder!
newfolder="chromopainterparsedout"

# call(["mkdir %s" %(newfolder)],shell=True)


# the new file is going the be created from the input file, and stored in the newfolder.
gm=[]
# gmfile=newfolder+"/"+name+'_parsedGM'
gmfile=newfolder+"/"+chrom+'_parsedGM'


print gmfile
gmfile=open(gmfile,'w')


for i in s:
        #print i
	if 'HAP' in i:
		split=i.replace('\n','').split(' ')

		hap=split[2]
	else:
		try:
			linecode=int(i[0])
			if int(i[0]) == 2: # chromopainter produces several samples of painted chromosomes. By default 10, in case only one sample is produced, take the first. There is a problem because the default is 10. Then 1 and 10 are confounded, use at least 2 samples
				gm.append(" ".join([hap,i[2:]]))
				# gm.append(" ".join([hap,i[2:50]])) # if you want a subset of all positions just to create a toy output, uncomment this and commend the upper line.
		except ValueError:
			pass   
for i in gm:
	#print i[0:10]
	gmfile.write(i)

gmfile.close()
