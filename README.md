# README.md

## Background

DBSLMMinterval is a software product for estimating polygenic scores in large GWAS data sets, such as those in the UK Biobank. We use the DBSLMM method for calculating point predictions for every subject's polygenic score. 

To construct prediction sets for quantitative traits, we first calculate the $$Var(\hat{\tilde y})$$ under the DBSLMM model. 

## Requirements

DBSLMMinterval requires installations of PLINK1 and PLINK2. It is also useful to have
GEMMA and LDscore installed locally. 

## Implementation

Unlike DBSLMM, DBSLMMinterval relies on a "control file" to specify inputs.

### Control file structure

The control file is a comma-separated values text file with one row per chromosome. 

Because the underlying code iterates over the chromosomes, we require that users 
enter in the control file eleven values per chromosome: 

1. Chromosome number  
2. b, the block information directory path  
3. eff, the effects output file path  
4. h, heritability estimate, a number between 0 and 1.    
5. l, summary file path for large effect markers  
6. mafMax, maximum discrepancy in MAF between reference and observed data, a number between 0 and 1      
7. n, sample size of the observed data, a positive integer    
8. nsnp, number of SNPs in the observed data, a positive integer  
9. r, file path for the reference panel bfile    
10. s, summary file path for small effect markers  
11. t, thread number, an integer, 1 or greater  

In the above, "observed data" refers to the data that is not the reference panel. 

Quotation marks must not be used in the control file.  

Control files must not have a header row. That is, the first row must be the information for Chr 1. 


