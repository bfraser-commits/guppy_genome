# Scripts for generating coverage data for Fraser et al. 2020
### Step 1 - Estimate coverage from bams
  * First step takes a directory of BAMs and calculates coverage using deeptools bamcoverage (https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html).
  * Coverage is calculated here in 50 BP windows anad averaged over 75 BP.
  * Outputs bedfile of coverage estimates per individual.
  * Script is split into full chr and scafs to account for different blacklisted regions.

### Step 2 - Mash individuals to population bedfiles
  * Takes as input a directory of individual bedfiles and calculates, per chromosome/scaffold, the population averaged coverage.
  * Populations here are males and females of each sampling population (6 populations, 3 rivers).
  * Includes additional filtering step that removes windows with coverage estimates 4 times greater than the normalised mean (1.0).
  * Outputs bedfile of coverage estimates per population per chromosome.
  
### Step 3 - Non-overlapping sliding window bedfiles
  * Takes as input directory of population-averaged chromosome coverage bedfiles and averages over the genome into user-defined non-overlapping windows.
  * Outputs final bedfiles.
  
### Step 4 - Comparison of male-female coverage
  * Rscript takes bedfiles from step 3 and inspects male/female coverage across the genome.
  * Includes production diagnostic figures used to evaluate per chromosome male/female coverage ratios.
  * Outputs male/female coverage comparisons for further plotting and analysis alongside population genomics data.
