# Minimac4ImputationAssessment
Imputation of GWAS data with a reference panel using Minimac4. An assessment of the imputation accuracy is also performed.

Input: target file in vcf format, reference file in vcf format
Output: Imputed target file, text file of imputation accuracy
Requirements: Python 3.5 or higher


Typical command line usage
python3 Minimac4ImputationAccuracy.py --tar target.vcf.gz --ref reference_panel.vcf.gz --chr 1 --cpus 8
Please separate your files by chromosome before running. At the moment there is no implementation to do this for you.

