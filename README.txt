########### Ancestry GWA #############
###### Moises Exposito-Alonso ########
### moisesexpositoalonso@gmail.com ###
####### Max Planck Institute #########
############ Weigel Lab ##############
 

## State the problem
Quantitative genetics is a well fundamented science that usually focuses on the effect of specific alleles (ACGT) in a quantitatie phenotype. 
Oftentimes achieving this task is challenging due to confounding structure in the data such as historical population conections and complex
admixture. The solution was to correct for it using a random variable with a kinship matrix in a mixed model framework. This has the advantage
of reducing false positive SNPs, but imposes the limitation of missing adaptive changes linked to population history processes, a concept of
high interest in evolutionary biology but of no relevance for quantative genetics in medicine or breeding. 

## State of the art
Use the concept of genome wide screening of GWAs and the concept of ancestry assignation from softwares such as Hapmix and Chromosome painter. 
Instead of associating allele states with a phenotype, we rather associate ancestry states at each position in the genome. 

## My analyses

What I have of 762 individuals: (1) a phenotype, drought resistance. (2) vcf, plink, etc of genomes. (3) Chromosome painter analyses. Unguided 
or all to all mode. (4) genetic clusters assigned by ADMIXTURE.

## Description of files


## Pipeline steps and files description

(A) Parse the output from Chromosome painter. There is one .sample file per chromosome with several painted chromosomes runs per indidivual. Provided 
the file name, this script will produce another, more readable, file named '_parsedGM' at the end. This has to be done for each chromosome. Make sure 
that the .samples.out file only has a number that refers to the chromosome, otherwise there will be problems downstream.

E.g.run
python parse_sampleout_file_otherfolder.py /ebio/abt6_projects8/ath_1001G_history/finestructure/guideddrought/chr1.samples.out


(B) Run a chromopainted gwa. Since I know where my results are, I just give as an argument the chromosome number. I need to provide a label of the input, 
which are the genomes.
This is the label input, is inside the R file, but I modify in great deal to do several analyses.
/abt6_projects9/ath_1001G_image_pheno/experiment_218_droughtgwa/fineancestry/label-input-relicts.csv
Also need to provide a phenotype. In this case I have my drought phenotype (and others) in the file
../MASTER_DATA.RData
E.g. run
Rscript finestructure_descriptivegwa_pipe.R 1
Rscript ancestrygwa.R

Inside this file, there is a key function, fineANOVA, which after using apply in sampleout, will merge the line of ancestry and the phenotype file and 
compute a p-value of the differences of phenotypes across ancestry haplotypes using a kruskal wallis test. If the phenotype is normally distributed 
also ANOVA could be done, but I implemented a non-parametric test in case the variance of the phenotype per groups is also difficult (non-homocedastic).

The output p-value list is saved as an object file named finegwa_chrX.RObject

This needs an extra label input file to associate the chromopainter output with phenotypes and population structure. The first column should be the
genome identifiers, the second the population group (discrete variable) or principal component as one possible descriptor of the structure (continuous
variable).

 
(C) Generate a genome map. From the haplotypes input of chromosomepainter, the python file extracts the row of positions, and the R file puts all the 
chromosome positions and generates an R object called chrpos.RObject.

E.g. run
python parsie_positions_fromhaplotypes.py 1
parse_positions_chr.R

(D) Posthot analyses
(D.1) Plot a classic Manhattan plot. Given a fingewa_chrX.RObject this file will load the finegwas results of the 5 genomes and the genome map positions and 
produe a manhattan plot.

E.g. run
finestructure_gwa_plot_ALL.R

(D.2) Get the top SNPs. This file will produce a top list of your SNPs. A given number of top hits per chromosome and write a table finegwa_TOP_xname.csv

E.g. run
finestr *** to be added

(D.3) Painted chromosome and the p-value from analyses in one plot. The scripts goes over the top SNP table and plots a given number of SNPs up and 
downstream from the SNP hit.

E.g. run
painted_chromosome_hits.R

(E) Other miscelaneous analyses. 




