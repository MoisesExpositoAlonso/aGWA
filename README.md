------------------------------------------------------------

# Ancestry Genome Wide Association (aGWA) #
#### Moises Exposito-Alonso #
##### moisesexpositoalonso@gmail.com #
##### moisesexpositoalonso.wordpress.com #
##### Weigel Lab Max Planck Institute #

------------------------------------------------------------
 
These set of scripts are free to use, modify and share, but it comes with no warranty. Emails with bugs or questions: moisesexpositoalonso@gmail.com

If you use this methodology, please cite as:
M Exposito-Alonso. aGWA: genome wide association analysis with ancestry information. DOI:dx.doi  (soon will be added)


### WISHLIST OF FUTURE IMPLEMENTATIONS
* ADD EXAMPLE TOY DATA!!!
* Extend the wrapper scripts to run autonomously chromosome painter analyses. (so far the user has to run it *a priori*)
* Re-structure the scripts in a workflow management system such as snakemake


### State the problem & state of the art
Quantitative genetics is a well fundamented science that usually focuses on the effect of specific alleles (ACGT) in a quantitatie phenotype. Oftentimes achieving this task is challenging due to confounding structure in the data such as historical population conections and complex admixture. The solution so far has been mainly account for genome background incorporating a random variable with a kinship matrix in a mixed model framework. This has the advantage of reducing false positive SNPs, but imposes the limitation that any relevant genetic variant linked to population history processes will be removed, being these alleles of particular interest when studying the genetic basis of local adaptation. These analyses can also be applied when susceptibility to a disease varies between ethnical background in humans; some Labs have developed similar ideas and tools and called it admixture mapping.

Use the concept of genome wide screening of GWAs and the concept of ancestry assignation from currently widely used softwares such as Chromosome painter to assign ancestry by chunks in the genome. Instead of performing association analyses with specific allele stats (e.g. A *vs* C), do so with ancestry states defined by chromosome painter.

### This repository has been created to share the analyses from Exposito-Alonso et al. 201X (*to be submitted*)

The analyses were based in 762 individuals from the plant _Arabdisopsis thaliana_ belonging to the [1001 genomes project](1001g.org). For these we had (1) a continuous phenotype. (2) whole genome sequences. (3) Chromosome painter analyses from the unguided form (all to all). (4) genetic clusters assigned by ADMIXTURE or other arbitrary assignment of population (e.g. geograhpic, which actually coincides in _A.thaliana_)


### Pipeline steps and files description
All required commands are in the file sh_ANALYSES-ANCESTRYGWA.sh. This bash script has to be edited with the location of your chromosome painter .samples.out file, and the input for chromosome painter, the haplotype files (from this we extract the locations of each SNP). 
Then, as arguments for the sh_ANALYSES-ANCESTRYGWA.sh script, the (1st) must be a .tsv file with the information per genome (labelinput), the (2nd) the type of analysis, "discrete" (implemented so far), and (3rd) the name of this run, e.g. "run_1".

An example of how to run it is found in the file "example_job.sh", and would be like below. I added the "time" command to know how long it takes to run, and redirect the output to a log file named as the analysis.

``` sh
time bash sh-ANALYSES_ANCESTRYGWA.sh label_input_genomes.tsv discrete run_1 > run_1.log
```

The labelimput file should be a tab separated file with first column the id of the genome, second column the population, third column the continuous phenotype. It would look like this:

``` sh
159	9	-0.000591739079573827
265	9	0.000371342214314237
350	9	NA
351	9	NA
```

#### Steps in the analyses:
1. Parse the output from Chromosome painter. There is one .sample file per chromosome with several painted chromosomes runs per indidivual. Provided the file name, this script will produce another, more readable, file named 'parsedGM' at the end, and will be stored in the folder 'chromopainterparsedout/'. This has to be done for each chromosome. Make sure that the .samples.out file only has a number that refers to the chromosome, otherwise there will be problems downstream. This has to be edited by the user in the sh_ANALYSES-ANCESTRYGWA.sh file. 

 ```sh
 # example parsing Chromosome painter output:
 chromo1=chr1-fs.samples.out
 python ../parse_sampleout_file_otherfolder.py $chromo1 1 
 ```

2. Generate the position and chromosome map from the input used in Chromosome painter (haplotype files). This will generate files as '1_chrpositions'.

 ```sh
 # example parsing haplotype files to get SNP positions:
 haplo1='1135-imp-1.haplotypes-drought'
 python ../parse_positions_fromhaplotypes.py $haplo1 1 
 ```

3. Join the SNP positions and get ready for later analyses. This will generate R objects as 'chr_1.RObject' and later 'chrpos.RObject' for all chromosomes.

 ```sh
 # Example loading SNP positions in chromosome 1
 Rscript../parse_chr_separated.R 1 
 # After each chromosome is read, then join those in a single object
 Rscript ../parsejoin_chr_positions.R 1 2 3 4 5
 ```

4. Run a chromopainted gwa. It requires (1) the labels file, (2) the chromosome number (knowing the number will look for the parsed genome from above), and (3) what type of analysis (currently discrete). The test is based in a Generalized Linear Model (GLM), therefore the script is easily generalizable to phenotypes with distributions other than normal. This analyses will generate a set of p-values along the genome, but for each position also will calculate an R2, measuring how much phenotype variance the ancestry of the SNP explains. It will be stored in the folder "results/" under the name "agwa_chr1.RObject" for p-values and "agwa_r2_chr1.RObject" for R2 values. 

 ```sh
 labelinput=label_input_genomes.tsv
 typeagwa="discrete"
 # Example aGWA in chromosome 1
 Rscript ../ancestrygwa.R $labelinput 1 $typeagwa 
 # After each chromosome is read, then join those in a single object
 Rscript ../parsejoin_chr_positions.R 1 2 3 4 5
 ```

5. To reduce the false positives inherent to multiple testing and population structure, I designed an empirical p-value distribution that will be used to correct the p-values from section 4. This permutation procedure consists in shuffling SNPs within individual at a distance enough to get out of a painted chromosome. The permutation of a focal position with another position within an individual will allow to maintain the overall ancestry of the individual while disrupting a specific block that might have been inherited in multiple individuals, which might be the source of phenotypic variance and what is aimed to identify. Therefore the p-values calculated only would reflect the effect of background ancestry and capture the stockasticity of p-values when performing many tests.

  1. Characterizing the chunks size of chromosome painter. To be able to perform a sucessful shuffling of SNPs we need to determine the distribution of painted chunks in the genome. To capture an average for all the genomes and different regions of the chromosome, a walking program runs 1 million steps per genome, counting how many steps it takes to change ancestry.  Every time it changes ancestry, i.e. a new block, it jumps randomly somewhere else in the genome to be able to cover all the genomic heterogeneity. This saves a file under "chromopainterparsedout/" folder named as the parsed genome matrix with the ending "BLOCKSIZES", and plots the block length distribution per genome inside "plots/" named as "chr1_block_length_distribution.pdf". The block length distribution of all chromosomes joined is used in the 5.2 section as the density probability distribution of the distance at which a SNP is shuffled.
   
     ```sh
     # Example for chromosome 1:
     python ../block_length_dist.py 1
     # To visualize the block lengths there is a plotting script. It stores the plots at "plots/" subfolder.
     Rscript ../plot_block_length_dist.R 1   
    ```
 
  2. aGWA under permuted datasets. Each chromosome has 50 replicates by default. Each replicate consists in 100 shuffled SNPs. This will produce many objects under "results/" named by the chromosome number and the replicate number, for example: "rep1_agwa_permuted_chr1.RObject". Afterwards, all p-values are joined in the object "results/agwa_permuted_all.RObject". Once joined, a distribution of all empirical p-values is plotted as "plots/agwa_permuted_all.RObject.hist.pdf".
   
     ```sh
     # The script already has a loop to paralelize many anayses replicates (50 by default). 
     # Example for chromosome 1 and replicate 1. 
     Rscript ../permutation_master.R $labelinput 1 1 $typeagwa 
     # Then all p-values from the permutation procedures are joined. 
     Rscript ../permutation_join.R
     ```
   
6. Relativize aGWA p-values. Each p-value in section 4 is relativized as the percentile in the distribution produced in 5.2. That is, if the p value is 0.005 and in the distribution, which has a total of 25000 empirical p-values,there are 1000 below 0.005, the corrected p-value is 0.04. All aGWA p-values with no empirical p-values below, are assigned a corrected p-value of 1/25000 = 0.00004. Relativized agwa p-values are stored as "results/agwa_chr1.RObject_empiricalpval.RObject". A corresponding plot is generated for the uncorrected and the corrected p-values per chromosome in files under "plots" ending in ".hist.pdf"
  
  ```sh
   # Example for chromosome 1. 
   Rscript ../relativizepval.R 1  
   ```

7. To plot the results, we need to join the computed p-values with the corresponding chromosome positions. 
  
  ```sh
   # Run the next. This script already parses the "result/" folder to get all agwa and agwa empiricalpval objects.
   Rscript ../combine_pvalue_positions.R
   ```
   
8. Plot. A classic manhattan is plotted. There are two versions, one is a multiple plot with one chromosome one row (flag "perchromosome"). The other version is all chromosomes concatenated (flag "cumulative"). I prefer "perchromosome"
  
  ```sh
   # Run the next. By default it plots the empirical p-value object, since they regular agwa has a really high false discovery rate.  
   Rscript ../manhattan.R perchromosome
   ```
  * This produces two plots, one with all the SNPs that can be many Mb size, "plots/aGWAmanhattan.pdf", and a subsample of 10000 points, "plots/subsampleplotaGWAmanhattan.pdf". 

9. Get top SNPs. A single script will retrieve all SNPs with p-values below a threshold (from 0 to 1), or a given number of top SNPs. A table is produced "tables/toptable_0.001.tsv", with the resulting SNPs. 

  ```sh
   # Run the next for all the SNPs with p-value<0.001 
   Rscript ../topsnps.R 0.001 
   # Run the next for the top 100 SNPs 
   Rscript ../topsnps.R 100 
   ```

### NOTES

Every R scripts sources a library, ancestrygwa_functions.R, that contains all the functions. 

Debugging should be fairly easy since steps are very well defined and each runned command is piped to store the standard output and standar error to a log file, all stored under the folder "logs/".



