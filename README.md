# Minimac4 Imputation Pipeline

In order to obtain higher quality GWAS data, imputation software can be used to fill in any missing snps. Of course, there is always the problem that the imputed values are incorrect and can be detrimental to a GWAS analysis. This pipeline will filter duplicate and incorrectly coded entries snps, prephase, impute, and analyze the output to determine the accuracy of the imputation. It is designed to work on at a whole genome level on chromosomes 1-22.

## Getting Starated

### Prerequisites

* Built and tested on Ubuntu 16.04.6 LTS with Bash terminal. There is no gurantee this can run on other systems.
* Python 3.5 or higher.
* VCF files for target and reference data, both GRCh38 build, bgzipped with a tabix index.
* The following python third party packages.

```python
pip install pandas
pip install numpy
pip install matplotlib
```

* The following program directories need to be on path.

```bash
htslib
vcflib
vcftools
Minimac3
Minimac4
```

### Installing

Download this repository with the following command.

```Git
git clone https://github.com/FischyM/Imputation-Pipeline.git
```

Decompress the genetic map files in the map directory.

```Bash
cd map
unzip plink.GRCh38.map.zip
gunzip genetic_map_hg38_withX.minimac.txt.gz
```

To see the usage statement, run the program with the '-h' flag.

```bash
python Minimac4ImputationAccuracy.py -h
```

```bash
usage: Minimac4ImputationAccuracy.py [-h] -t TAR -r REF [--cpus CPUS]

Minimac4ImputationAccuracy.py -t <target file> -r <reference file> -cpu
<default=4>

optional arguments:
  -h, --help         show this help message and exit
  -t TAR, --tar TAR  target GWAS dataset as vcf.gz
  -r REF, --ref REF  reference dataset as vcf.gz
  --cpus CPUS        number of threads to use

```

## Authors

Mathew Fischbach, Undergraduate Researcher of Myers Lab: <http://csbio.cs.umn.edu/people.html>

## Acknowledgements

Pipeline implementation and workflow based on the work of Xiaotong Liu, Myers Lab.
