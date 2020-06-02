# Minimac4 Imputation Assessment
Imputation of GWAS data with a reference panel using Minimac4. An assessment of the imputation accuracy is also performed.

Input: target file in vcf format, reference file in vcf format
Output: Imputed target file, text file of imputation accuracy
Requirements: Python 3.5 or higher. Please have the following program executables in your PATH environment: htslib (includes tabix, bgzip), vcflib, vcftools, Minimac3, Minimac4. Pre-phasing using beagle is within this repository as well as genetic map files found in the map directory. They will need to be decompressed before proceeding. 


Typical command line usage
python3 Minimac4ImputationAccuracy.py --tar target.vcf.gz --ref reference_panel.vcf.gz --cpus 8
This program will separate your files by chromosomes using bgzip and tabix and continue with the imputation and analysis one chromosome at a time. Results of the imputation and its analysis can be found in the analysis directory under each specified chromosome. 


