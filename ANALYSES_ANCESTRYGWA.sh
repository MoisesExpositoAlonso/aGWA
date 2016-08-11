#!/bin/bash


## This bash script has the list and order of commands that need to be run.
# python parse_sampleout_file_otherfolder.py finestructure/weiout/chr1-fs.samples.out
# python parse_sampleout_file_otherfolder.py finestructure/weiout/chr2-fs.samples.out
# python parse_sampleout_file_otherfolder.py finestructure/weiout/chr3-fs.samples.out
# python parse_sampleout_file_otherfolder.py finestructure/weiout/chr4-fs.samples.out
# python parse_sampleout_file_otherfolder.py finestructure/weiout/chr5-fs.samples.out

## Create the SNP map
# python parse_positions_fromhaplotypes.py 1135-imp-2.haplotypes-drought
# python parse_positions_fromhaplotypes.py 1135-imp-3.haplotypes-drought
# python parse_positions_fromhaplotypes.py 1135-imp-4.haplotypes-drought
# python parse_positions_fromhaplotypes.py 1135-imp-1.haplotypes-drought
# python parse_positions_fromhaplotypes.py 1135-imp-5.haplotypes-drought

# parse_positions_chr.R

## Now the ancestry GWA analysis
# Rscript ancestrygwa.R label_input_genomes.tsv 1 & 
# Rscript ancestrygwa.R label_input_genomes.tsv 2 &
# Rscript ancestrygwa.R label_input_genomes.tsv 3 &
# Rscript ancestrygwa.R label_input_genomes.tsv 4 &
# Rscript ancestrygwa.R label_input_genomes.tsv 5 

# for c in 1 2 3 4 5
# do 
# namelog="agwa_"$c".log"
# echo $namelog
# Rscript ancestrygwa.R label_input_genomes.tsv $c & > $namelog &
# done
# wait

## Get the average length of blocks

# python block_length_dist.py chromopainterparsedout/chr1-fs.samples.out_parsedGM
# python block_length_dist.py chromopainterparsedout/chr2-fs.samples.out_parsedGM
# python block_length_dist.py chromopainterparsedout/chr3-fs.samples.out_parsedGM
# python block_length_dist.py chromopainterparsedout/chr4-fs.samples.out_parsedGM
# python block_length_dist.py chromopainterparsedout/chr5-fs.samples.out_parsedGM

# plot_block_length_dist.R

## Generate empirical distribution

# Rscript empirical_distribution_bypermutation.R label_input_genomes.tsv 1 &
# Rscript empirical_distribution_bypermutation.R label_input_genomes.tsv 2 &
# Rscript empirical_distribution_bypermutation.R label_input_genomes.tsv 3 &
# Rscript empirical_distribution_bypermutation.R label_input_genomes.tsv 4 &
# Rscript empirical_distribution_bypermutation.R label_input_genomes.tsv 5 

echo'Generate empirical distribution'
for c in 1 2 3 4 5
do 
namelog= "permut_"$c".txt"
echo $namelog
Rscript empirical_distribution_bypermutation.R label_input_genomes.tsv $c & > $namelog &
done
wait

# Rscript relativizepval.R 1 &
# Rscript relativizepval.R 2 &
# Rscript relativizepval.R 3 &
# Rscript relativizepval.R 4 &
# Rscript relativizepval.R 5 &

echo'relativize pvalues'
for c in 1 2 3 4 5
do 
namelog= "relpval_"$c".txt"
echo $namelog
Rscript relativizepval.R $c & > $namelog &
done
wait

## Combine output with a chromosome and position map

Rscript combine_pvalue_position.R  > combinepvalues.log ######## NEED TO ADAPT FROM RELATIVIZED P VALUE ALSO!



# run the gwa plot and top SNP

# Rstudio aGWA_manhattan.R

# top SNPs

# Rstudio topsnps.R

# Painted chromosome around hit

# Rstudio paintedhits.R

