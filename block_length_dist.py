import sys
import random

# get the input file, the painted chromosome from an argument
try:
 chrom=sys.argv[1]
 filename=chrom+"_parsedGM"
 filename="chromopainterparsedout/"+filename
except IndexError:
 print "no input file provided!"

print filename
painted=open(filename,"r")


# create a blocksizes file, where the different lengths are stored
newfile=filename+"_BLOCKSIZES"
blocksizes=open(newfile,"w")

# several counters for the chain
maxacccounter=762
maxposcounter=1000000

genomecounter=0
poscounter=0
breakscounter=0
realcounter=0


for row in painted:
 #split positions 
 row=row.split(' ')
 print "genome id number", row[0]

 #count how many rows (genomes) have been walked
 genomecounter=genomecounter+1
 print "genomes walked", genomecounter
 # break if already walked many genomes
 if genomecounter  >= maxacccounter:
  break
 
 #random starting point in the chromosome
 #startpos=2
 startpos=random.randrange(2, len(row) ) 

 realcounter=0
 for snppos in range(startpos,len(row) ): 
  #count how many iterations within a genome
  realcounter=realcounter+1
  
  # break if already walked many positions
  if realcounter >= maxposcounter:
#   print "positions walked", realcounter
   break

  # condition test: equal or different haplotype?
  if row[snppos] == row[snppos-1]:
#   print "equal?", row[snppos], row[snppos-1]
   poscounter=poscounter+1
  else:
#   print "different?", row[snppos], row[snppos-1]
   print "steps walked in a single haploblock ",poscounter
   blocksizes.write("%s\n" %(str(poscounter)))
   poscounter=0

blocksizes.close()



