# About
This repository contains a collection of scripts used in the "Ancient human genomes reveal a hidden history of strong selection in Eurasia" Paper. For additional information about using the scripts please contact:

Yassine Souilmi: yassine.souilmi@adelaide.edu.au
Raymond Tobler: raymond.tobler@adelaide.edu.au


## Notes on data format
The bash wrapper script that runs creates SweepFinder2 input (SFS) expects the data in Plink format. You can get that by using convertf utility available as part of the AdmixTools suite. However, make sure that the .fam is properly formated with the population IDs in the first column, and that the phenotype column is setup to either 1, or 0 and not -9. Also, sex chromososmes are ignored and can be removed from the dataset.

Also, needed is a refernce file (fasta). Plink has the upleasant habit of polarising each SNP based on the minor allele frequency whithin the dataset. To remedy that prior to the generation of the SFS, we repolarise each SNP using the reference allele using bcftools (load the module if necessary).

In the shared repository there's an example file for a subset of the 1240k SNP set. As the polarisation step could be accomplished with Plink 2.x if used (--update-ref-allele refAllele.txt).

## Running SweepFinder2
The warapper (sweepfinder2_splitByPop.sh under scripts folder)script is a simple bash that creates an SFS file for each chromosome for each population in the plink file. Also, make sure to edit the file to point to the righ location of the vcf2SF.py python script (scripts folder), which depend on python3 The script could be used as follows:

```
$ bash sweepfinder2_splitByPop.sh <plink file basename> <ref.fasta>
```

**Note:** Please note that sweepfinder jobs (per chromosome, per population) takes a fair amount of time to run, and scales up principally with the number of loci. Each job takes for 1M SNPs between 8-12h walltime. So it is advised to replace the SweepFinder2 calls (in the bash script) with a job submission. And/or replace the loop inside the script with parallel or xargs to activate parallelisation. Each job requires 2 cores and 2GB of ram.

**Note2:** For example input and output files please check: [https://github.com/ACAD-UofA/PolyLink/tree/master/Example](https://github.com/ACAD-UofA/PolyLink/tree/master/Example)

## R scripts
- `Determine_outlier_genes_and_sweeps.R`: helps assing SweepFinder2 scores (CLR) and determine outlier sweeps for each population.
- `Fst_analyses.R`: helps compute *Fst* by implementing the Weir-Cockerham estimator.
- `OutFLANK.R`: fast implementation of the OutFlank method to determin *Fst* outliers.
