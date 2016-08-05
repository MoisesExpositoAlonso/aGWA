#!/bin/bash


# This bash script has the list and order of commands that need to be run.

python parse_sampleout_file_otherfolder.py chr1.samples.out
python parse_sampleout_file_otherfolder.py chr2.samples.out
python parse_sampleout_file_otherfolder.py chr3.samples.out
python parse_sampleout_file_otherfolder.py chr4.samples.out
python parse_sampleout_file_otherfolder.py chr5.samples.out

# Create the SNP map
python parsie_positions_fromhaplotypes.py 1
python parsie_positions_fromhaplotypes.py 2
python parsie_positions_fromhaplotypes.py 3
python parsie_positions_fromhaplotypes.py 4
python parsie_positions_fromhaplotypes.py 5

parse_positions_chr.R


# Not the ancestry GWA analysis

Rscript ancestrygwa.R label_input_genomes.tsv 1 
Rscript ancestrygwa.R label_input_genomes.tsv 2
Rscript ancestrygwa.R label_input_genomes.tsv 3
Rscript ancestrygwa.R label_input_genomes.tsv 4
Rscript ancestrygwa.R label_input_genomes.tsv 5

# Combine output with a chromosome and position map

Rscript combine_pvalue_position.R


# echo ' prepare label files ' 
# Rscript adapt_files_finestructure_onlyasians.R
# Rscript adapt_files_finestructure_onlyrelicts.R
# Rscript adapt_files_finestructure_onlyswedes.R
# Rscript adapt_files_finestructure_onlymediterraneans.R
# # this runs chromosomepainting
# cd finestructure
# echo 'running guided painting for drought genotypes ...'
# bash guidedpainting_pipe.sh onlysweds &
# bash guidedpainting_pipe.sh onlyrelicts &
# bash guidedpainting_pipe.sh onlyasians &
# bash guidedpainting_pipe.sh onlymediterraneans &
#bash guidedpainting_drough.sh

# cd ..
# # this prepares the output from chromosomepainting and runs an ANOVA GWA.
echo 'run finegwa analyses on each painted chromosome ...'
# bash run-finegwa_ONLYPARSE.sh
#bash run-finegwa.sh > run-finegwa.log

# echo ' get the positions of chromosomes'
# # preapre positions for the gwa plot and top SNPs
# Rscript parse_positions_chr.R

# # run the gwa plot and top SNP
echo ' plot finegwa ...'
Rscript finestructure_gwa_plot_ALL.R res > finestructureplotall.log &
Rscript finestructure_gwa_plot_ALL.R  resNSwed > finestructureplotall_resNSwed.log &
Rscript finestructure_gwa_plot_ALL.R resMed > finestructureplotall_resMed.log &
Rscript finestructure_gwa_plot_ALL.R resRel > finestructureplotall_resRel.log &
Rscript finestructure_gwa_plot_ALL.R resCat > finestructureplotall_resCat.log &

# # get the first top 5 snps all chromosomes
# echo 'plot top hits painted chromosomes'
#Rscript painted_chromosome_hits.R > paintedhits.log



