# Impute5 Imputation Pipeline

In order to obtain higher quality GWAS data, imputation software can be used to fill in any missing snps. Of course, there is always the problem that the imputed values are incorrect and can be detrimental to a GWAS analysis. This pipeline will filter duplicate and incorrectly coded entries snps, prephase, impute, and analyze the output to determine the accuracy of the imputation. It is designed to work on at a whole genome level on chromosomes 1-22.

## Getting Starated

### Prerequisites

* Built and tested on Ubuntu 16.04.6 LTS with Bash terminal. There is no gurantee this can run on other systems.
* Python 3.8
* VCF files for target data, GRCh38 build, bgzipped with a tabix index.
* See requirements.txt for python modules
* Required map files for Impute5 and Shapeit4, provided here
* The following program directories need to be on path.

```bash
htslib
bcftools
impute_v5.1
shapeit4
```

### Installing

Download this repository with the following command.

```Git
git clone https://github.com/FischyM/Imputation-Pipeline.git
```

Decompress the genetic map files in the map directory.

```Bash
cd map/
tar -xvzf genetic_maps.b38.tar.gz
```

To see the usage statement, run the program with the '-h' flag.

```bash
./Impute5Pipeline.py -h

usage: Impute5Pipeline.py [-h] -t TARHAPS -r REFHAPS [-p PREFIX] [-c CPUS] [-m MEMORY]

optional arguments:
  -h, --help            show this help message and exit
  -t TARHAPS, --tarHaps TARHAPS
                        Target GWAS dataset in vcfn vcf bgziped, or bcf format as a singular file. This program will do the
                        separation
  -r REFHAPS, --refHaps REFHAPS
                        File list containing reference datasets (in order of chromosomes 1-22, one per line) in imp5 format.
  -p PREFIX, --prefix PREFIX
                        Prefix for output files.
  -c CPUS, --cpus CPUS  Number of threads to use. Suggested around 12 cores with approx 100 Gb RAM.
  -m MEMORY, --memory MEMORY
                        length of region to impute at a time. 10 cM is default. To reduce memory usage, use less than 10000000.

```

By default the shebang line is:
```python
#!/usr/bin/env python3.8
```

## Authors

Mathew Fischbach, Undergraduate Researcher of Myers Lab: <http://csbio.cs.umn.edu/people.html>

## Acknowledgements

Pipeline implementation and workflow based on the work of Xiaotong Liu, Myers Lab.
Thanks to Wen Wang for her guidance and help.
Thanks to Chad Myers and the rest of the Myers lab for their help and resources.
