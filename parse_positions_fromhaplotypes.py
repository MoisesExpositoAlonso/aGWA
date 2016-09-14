import sys

# receive input file
try:
 hapfile=sys.argv[1]
except IndexError:
 print "no input file"
 print "run as: python parse_positions_fromhaplotypes.py 1135-imp-1.haplotypes-drought 1" 
 exit()
hap=open(hapfile,'r')

# receive chrom number
try:
 chrom=sys.argv[2]
except IndexError:
 print "no input chromosome name"

hapout =open(''.join([chrom,"_chrpositions"]),'w')

counter=0
for line in hap:
    print line[0]
    counter =counter+1 
    
    if counter== 3:
        
    #if line[0] == 'P':
	hapout.write(line[2:])
        hapout.close()
        #print line
        print line [0:5]
        print line[2:5]
    elif counter > 3:
	break
    else:
        pass

print('done !')
