
# example run; python parse_positions_fromhaplotypes.py 1135-imp-1.haplotypes-drought   
import sys

# hapfile=open('toyhap.haplotypes','r')
hapfile=sys.argv[1]
hap=open(hapfile,'r')
#print hap
hapout =open(''.join([hapfile,"_chrpositions"]),'w')

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
