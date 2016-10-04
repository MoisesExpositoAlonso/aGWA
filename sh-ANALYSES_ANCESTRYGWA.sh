#!/bin/bash
echo 'interpreted arguments:'
echo $0
#============================================================================================#
## Define general variables ##

# Sample out (provide absolute paths)
chromo1='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/weiout/chr1-fs.samples.out'
chromo2='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/weiout/chr2-fs.samples.out'
chromo3='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/weiout/chr3-fs.samples.out'
chromo4='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/weiout/chr4-fs.samples.out'
chromo5='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/weiout/chr5-fs.samples.out'

# Haplotype files (provide absolute paths)
haplo1='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/1135-imp-1.haplotypes-drought'
haplo2='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/1135-imp-2.haplotypes-drought'
haplo3='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/1135-imp-3.haplotypes-drought'
haplo4='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/1135-imp-4.haplotypes-drought'
haplo5='/ebio/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/finestructure/1135-imp-5.haplotypes-drought'


# Label input of structure and phenotypes (no absolute path, ony name)
labelinput=$1
#labelinput='label_input_genomes.tsv'
labelinput='../'$labelinput

# The type of aGWA you want
#typeagwa='continuous' # uncomment to use
#typeagwa='discrete' # uncomment to use
typeagwa=$2

## Iterations of empirical p value computation. The following number of iterations are done per chromosome. Each iteration is
# comprised by 100 tests of 100 shuffled SNP positions. Careful with high number of interations since they are done in parallel
# and will start many processes in the system simultaneously.
iterations=50

# output folder name
#out='trialcontinuous'
out=$3

mkdir $out
cd $out ### MOVING THE WORKIND DIRECTORY !

## Generate some required folders
mkdir plots
mkdir logs
mkdir chromopainterparsedout
mkdir results
mkdir tables

# chromosome list (if you have many chromosomes, maybe you want to generate an array object within the bash script, like this. I don't use this object though)
chrlist=(1 2 3 4 5)
chr=${chrlist[@]}


#============================================================================================#
## Run analyses ##
# This can be runned separatedly. You can skip steps or run separatedly as long as you provide with the correct input.

## This bash script has the list and order of commands that need to be run.
echo "parse sampleout ..."
python ../parse_sampleout_file_otherfolder.py $chromo1 1 &
python ../parse_sampleout_file_otherfolder.py $chromo2 2 &
python ../parse_sampleout_file_otherfolder.py $chromo3 3 &
python ../parse_sampleout_file_otherfolder.py $chromo4 4 &
python ../parse_sampleout_file_otherfolder.py $chromo5 5 &
echo "... done parse sampleout"

## Create the SNP map
echo "parse positions of each SNP in chromosome..."
python ../parse_positions_fromhaplotypes.py $haplo1 1 &
python ../parse_positions_fromhaplotypes.py $haplo2 2 &
python ../parse_positions_fromhaplotypes.py $haplo3 3 &
python ../parse_positions_fromhaplotypes.py $haplo4 4 &
python ../parse_positions_fromhaplotypes.py $haplo5 5 &
echo "done parse positions"

wait

echo 'joining positions SNPs...'
for c in 1 2 3 4 5
do 
Rscript ../parse_chr_separated.R $c &
done
wait
#Rscript ../parse_positions_chr.R 1 2 3 4 5
Rscript ../parsejoin_chr_positions.R 1 2 3 4 5
wait
echo 'finished joined SNP map'

# Now the ancestry GWA analysis
echo "ancestry GWA scripts ...."
for c in 1 2 3 4 5
do 
namelog="logs/agwa_"$c".log"
echo $namelog
Rscript ../ancestrygwa.R $labelinput $c $typeagwa &> $namelog &
done
wait
'...finished aGWA scripts'


## Get the average length of blocks
echo 'calculate block length distribution'
for c in 1 2 3 4 5
do 
  namelog="logs/blocklength_"$c".txt"
  python ../block_length_dist.py $c &> $namelog &
done
wait

Rscript ../plot_block_length_dist.R 1 2 3 4 5 &
echo '...finished calculation of block length distribution'

## Generate empirical distribution
echo 'Generate empirical distribution...'
for c in 1 2 3 4 5
do 
	for r in  $(seq $iterations) # modify number of iterations by changing the variable
		do
		namelog="logs/permut_"$c"_rep_"$r".txt"
		echo $namelog
		Rscript ../permutation_master.R $labelinput $c $r $typeagwa &> $namelog & 
		done
  wait # uncomment this if the paralelization is too much!
	done
wait

echo '...joining empirical distributions...'
Rscript ../permutation_join.R  ### uncomment for genome wide empirical pvalue correction
echo ' ... finished empirical distribution calculation'
# echo '...joining empirical distributions by chromosome...'
# Rscript ../permutation_join_bychr.R 1 2 3 4 5
# echo ' ... finished empirical distribution calculation'


# relativize p values!
echo 'relativize pvalues ... '
for c in 1 2 3 4 5
do 
namelog="logs/relpval_"$c".txt"
echo $namelog
Rscript ../relativizepval.R $c &> $namelog &   ### uncomment for genome wide empirical pvalue correction
# Rscript ../relativizepval_bychr.R $c &> $namelog &
done
wait
echo '... finished relativize pvalues'


## Combine output with a chromosome and position map
echo 'combine p values and positions...'
Rscript ../combine_pvalue_position.R  > combinepvalues.log 
echo '...finished combination p values and positions'

# run the gwa plot and top SNP

Rscript ../manhattan.R perchromosome

## top SNPs
echo 'create top SNPs table ... threshold 0.001 by default'
Rscript ../topsnps.R 0.001


## Get annotations
# Rscript ../getannotations.R table/toptable_0.001.tsv

## Painted chromosome around hit
## ... you probably have to tune this one a bit.
#echo 'create plots of painted hits from results table'
#Rscript ../paintedhits.R $labelinput $typeagwa 200 toptable_0.001.tsv


