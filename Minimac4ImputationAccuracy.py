#!/usr/bin/python
# Copyright Mathew Fischbach, Myers Lab http://csbio.cs.umn.edu/, 2020
# code based on the work of Xiaotong Liu from Myers Lab

# Purpose:
# This program will take plink files and impute them with Minimac4.
# Analysis of the imputation will be done by minimac_info_analysis.py
# This program imputes on a single chromosome at a time. A simple bash script can
# be created to impute chromosomes 1 through 22

# Requirements: 
# Input: vcf file containing target data and vcf reference panel

import sys
import argparse
import os
import subprocess

startDir = os.getcwd()

def setUp(argv, targetFileName, refFileName, cpus):
    print('\nChecking input')

    parser = argparse.ArgumentParser(description='Minimac4ImputationAccuracy.py \
        -t <target file> -r <ref file> -chr <chromosome> -cpu <threads>\n')
    parser.add_argument('-t', '--tar', required=True, help='target GWAS dataset as vcf.gz')
    parser.add_argument('-r','--ref', required=True, help='reference dataset as vcf.gz')
    #parser.add_argument('-c', '--chr', default='None', help='chromosome to analyze')
    parser.add_argument('--cpus', '--cpu', dest='threads', type=str, default='4', help='number of threads to use')
    args = parser.parse_args()
    
    print('Successful input\n')
    targetFileName = args.tar
    refFileName = args.ref
    # chromosome = args.chr
    cpus = args.threads
    return targetFileName, refFileName, cpus

def runCommandLine(commands):
    try:
        process = subprocess.Popen(commands, shell=True)
        process.communicate()
    except subprocess.CalledProcessError as err:
        print('Error: ', err)
        sys.exit(1)

def checkPrograms():
    print('Checking required programs and directories')
    # set all programs used to the env PATH variable, except for beagle pre-phasing
    if not os.path.isdir('analysis'):
        print('directory "analysis" not found. Creating directory')
        os.mkdir('analysis')
        for i in range(1, 23):
            os.mkdir('analysis/chr' + str(i))
        
    if not os.path.isdir('beforePhasing'):
        print('directory "beforePhasing" not found. Creating directory')
        os.mkdir('beforePhasing')
        os.mkdir('beforePhasing/tar')
        os.mkdir('beforePhasing/ref')
        
    if not os.path.isdir('afterPhasing'):
        print('directory "afterPhasing" not found. Creating directory')
        os.mkdir('afterPhasing')
        os.mkdir('afterPhasing/tar')
        os.mkdir('afterPhasing/ref')

    if not os.path.isdir('map'):
        print('directory "map" not found, be sure to download from our GitHub and unzip plink.GRCh23.map.zip')
        sys.exit(1)
        
    if not os.path.isfile('beagle.27Apr20.b81.jar'):
        print('Cannot find beagle.27Apr20.b81.jar. Please download it to the current working directory')
        sys.exit(1)

    if 'htslib' not in os.environ.get('PATH'):
        print('Cannot find htslib on PATH')
        sys.exit(1)

    if 'vcftools' not in os.environ.get('PATH'):
        print('Cannot find vcftools on PATH')
        sys.exit(1)

    if 'Minimac3' not in os.environ.get('PATH'):
        print('Cannot find Minimac3 on PATH')
        sys.exit(1)

    if 'Minimac4' not in os.environ.get('PATH'):
        print('Cannot find Minimac4 on PATH')
        sys.exit(1)

    if 'vcflib' not in os.environ.get('PATH'):
        print('Cannot find vcflib on PATH')
        sys.exit(1)

    print('Successful, all programs/directories are present\n')

def separateFile(file, type):
    print('beginning to separate {f}'.format(f=file))
    if file[-2:] == 'gz':
        command = ['htsfile ' + file]
        try:
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
            stdOut = process.communicate()
        except subprocess.CalledProcessError as err:
            print('Error: ', err)
            sys.exit(1)
        if 'BGZF-compressed' not in stdOut[0].decode("utf-8"):
            print('{f} was not compressed by bgzip'.format(f=file))
            print('In order to efficiently divide this file, it needs to be decompressed and then recompressed with bgzip')
            print('decompressing now')
            runCommandLine(['gunzip ' + file]) # 
            filegz = file
            file = file[:-3] # remove .gz
            print('compressing vcf file with bgzip')
            runCommandLine(['bgzip ' + file])
        else:
            print('{f} is bgzipped already'.format(f=file))
        print('forcing a new indexing of vcf.gz file with tabix')
        filegz = file # to keep consistent
        runCommandLine(['tabix -f -p vcf ' + filegz])
    elif file[-3:] == 'vcf':
        filegz = file + '.gz'
        print(file)
        print('compressing vcf file with bgzip')
        runCommandLine(['bgzip ' + file])
        print('indexing vcf.gz file with tabix')
        runCommandLine(['tabix -f -p vcf ' + filegz])
    
    os.chdir(startDir + '/beforePhasing/' + type)
    print('separating chromosomes')
    for i in range(1,23):
        chr = str(i)
        print('separating into chromosome {c}'.format(c=chr))
        runCommandLine(['tabix ../../' + filegz + ' -h ' + chr + ' > chr' + chr + '.vcf'])
    os.chdir(startDir)

def removeDuplicates(chr): # remove duplicate entries. Necessary for pre-phasing and should be in vcf format
    print('Beginning to remove duplicate entries')
    runCommandLine(['vcfuniq beforePhasing/tar/chr' + chr + '.vcf > \
        beforePhasing/tar/chr' + chr + '.noDups.vcf'])
    print('Done, duplicate entries removed\n')
    # here I could gzip it for later
    # out: beforePhasing/tar/chr#.noDups.vcf

def removeBadAllelesWithAWK(chr):
    print('Creating new target file removed of bad alleles')
    toCommandLine = ["awk '$1 ~ /^#/ {print $0;next} {if ($4 ~ /A|C|T|G/ && $5 ~ /.|A|C|T|G/) print $0}' \
    beforePhasing/tar/chr" + chr + ".noDups.vcf > beforePhasing/tar/chr" + chr + ".prephaseReady.vcf"]
    runCommandLine(toCommandLine)

def prePhase(chr, cpu):
    print('beginning prephasing with beagle')
    toCommandLine = ['java -jar beagle.27Apr20.b81.jar gt=beforePhasing/tar/chr' + chr + '.prephaseReady.vcf \
        out=afterPhasing/tar/prePhasedChr' + chr + \
        ' iterations=20 map=map/plink.chr' + chr + '.GRCh38.map \
        chrom=' + chr + \
        ' nthreads=' + cpu]
    runCommandLine(toCommandLine)
    print('Done prephasing')
    # outputs target as .vcf.gz

def convertRef(chr, cpu):
    # to run minimac4, ref file must be in m3vcf format. uses minimac3 to do this.
    # It was determined that we will let Minimac3 do its parameter estimation so that a map
    # file input is no longer needed. This makes the process easier when dealing with
    # different genome builds
    print('\nBeginning reference file conversion to m3vcf')
    toCommandLine = ['Minimac3 --refHaps beforePhasing/ref/chr' + chr + '.recode.vcf \
        --processReference \
        --prefix afterPhasing/ref/chr' + chr + ' --cpus ' + cpu + ' --rounds 5 --chr ' + chr]
    runCommandLine(toCommandLine)
    print('Done converting ref file to minimac3 format')
    # outputs as .m3vcf.gz

def impute(chr, cpu):
    # shouldn't need --ignore-duplicates since that was taken care of
    toCommandLine = ['minimac4 \
        --haps afterPhasing/tar/prePhasedChr' + chr + '.vcf.gz \
        --refHaps afterPhasing/ref/chr' + chr + '.m3vcf.gz \
        --prefix afterPhasing/imputed_chr' + chr + '\
        --referenceEstimates OFF \
        --mapFile map/genetic_map_hg38_withX.minimac.txt \
        --cpus ' + cpu]
    runCommandLine(toCommandLine)
    print('imputation using minimac4 - done')

def compare(chr):
    os.chdir('analysis/chr' + chr)
    runCommandLine(['python3 ../../minimac_info_analysis.py ../../afterPhasing/imputed_chr' + chr + '.info \
        > imputation_assessment.txt'])
    os.chdir(startDir)
 
    # simply sends stdout to text file 



######################################################################################
def main():
    # Set up argument inputs
    targetFileName = ''
    refFileName = ''
    cpus = 4
    targetFileName, refFileName, cpus = setUp(sys.argv[1:], targetFileName, refFileName, cpus)

    checkPrograms()
    separateFile(targetFileName, 'tar')
    separateFile(refFileName, 'ref')
    for i in range(1, 23):
        chromosome = str(i)
        removeDuplicates(chromosome)
        removeBadAllelesWithAWK(chromosome)
        prePhase(chromosome, cpus)
        convertRef(chromosome, cpus)
        impute(chromosome, cpus)
        compare(chromosome)

if __name__ == '__main__':
    main()

# current testing input
# python3 Minimac4ImputationAccuracy.py -t PD_NGRC_phs000196.GRCh38.chr1.vcf.gz -r ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.ver4.2.vcf.gz --cpus 10

### Unused functions ###
# def get1000Gen():
    # if not os.path.isdir('1000Genome'):
    #     print('1000Genome directory not, present. Creating directory')
    #     os.mkdir('1000Genome')
    #     os.chdir(startDir + '/1000Genome')
    #     for chr in range(1,23):
    #         runCommandLine(['wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr'\
    #             + str(chr) + '.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz'])
    #     os.chdir(startDir)
    #     if not os.path.isdir('map'):
    #         os.mkdir('map')
    #         os.chdir(startDir + '/map')
    #         print('getting GRCh38.map files')
    #         runCommandLine(['wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip'])
    #         runCommandLine(['unzip plink.GRCh38.map.zip'])
    #     os.chdir(startDir)

#3.2 No need to mask because we are using Minimac4 with it's .info file. 
# This function would need to be updated to work with an imputation program that does 
# not have the capabilities of Minimac4
# def maskTargetFile(chr):
    # runCommandLine(['vcftools/bin/vcftools --vcf afterPhasing/tar/chr' + chr + '.vcf \
    #     --freq --out afterPhasing/chr' + chr + '_to_mask'])

    # position_list = []
    # total = 0
    # chromosome = 0
    # with open('afterPhasing/chr' + chr + '_to_mask.frq') as f:
    #     read_data = f.readlines()
    #     for line in read_data: 
    #         words = line.split('\t')
    #         total += 1
    #         chromosome = words[0]
    #         position_list.append(words[1])

    # random_list = position_list
    # deleted_list = []
    # masked_number = 1000  # can also make it a percentage of total

    # i = 0
    # while i < masked_number:
    #     position = random.choice(random_list)
    #     deleted_list.append(position)
    #     random_list.remove(position)
    #     i += 1
    # deleted_list.sort()
    # print(deleted_list[0])

    # print('\ntotal number of snps: ' + str(total))
    # print('number of masked snps: ' + str(len(deleted_list)))
    # print('number of snps after masking: ' + str(len(random_list)) + '\n')

    # with open('afterPhasing/snps_after_masking.txt', 'w') as snps_after_mask:
    #     for j in random_list:
    #         temp_string1 = chromosome + '\t' + j + '\n'
    #         snps_after_mask.write(temp_string1)

    # with open('afterPhasing/snps_that_are_masked.txt', 'w') as snps_that_are_masked:
    #     for k in deleted_list:
    #         temp_string2 = chromosome + '\t' + k + '\n'
    #         snps_that_are_masked.write(temp_string2)

    # print("\ngathering of snps to create a randomly masked file and real file - done")

    # toCommandLine = ['vcftools/bin/vcftools --gzvcf ' + targetFileName + \
    #     ' --positions afterPhasing/snps_that_are_masked.txt --recode \
    #         --out afterPhasing/masked_snps_to_compare']
    # runCommandLine(toCommandLine)

    # runCommandLine(['gzip afterPhasing/masked_snps_to_compare.recode.vcf'])
    # print('\nreal snp mask file to compare - done')

    # toCommandLine = ['vcftools/bin/vcftools --gzvcf ' + targetFileName + \
    #     ' --positions afterPhasing/snps_after_masking.txt --recode \
    #         --out afterPhasing/snps_to_be_imputed']
    # runCommandLine(toCommandLine)

    # runCommandLine(['gzip afterPhasing/snps_to_be_imputed.recode.vcf'])
    # print('\nmasked snps to impute - done')

# def removeBadAlleles(chr): 
    # # takes snps from target and ref file and creates a new ref file from 1000Genome
    # print('finding list of snps from given target and reference files')
    # toCommandLine = ['vcftools \
    #     --vcf beforePhasing/tar/chr' + chr + '.noDups.vcf \
    #     --diff beforePhasing/tar/chr' + chr + '.noDups.vcf \
    #     --diff-site --out beforePhasing/snpPositionChr' + chr]
    # runCommandLine(toCommandLine) 

    # # create list of snps while removing poorly coded alleles 
    # # this is based off of personal experience. This can easily be changed
    # print('\nRemoving wrongly coded alleles')
    # snp_list = []
    # count = 0
    # with open('beforePhasing/snpPositionChr' + chr + '.diff.sites_in_files') as f:
    #     read_data = f.readlines()
    #     for i in range(1, len(read_data)): # skip header line
    #         words = read_data[i].split('\t')
    #         # Change the logic to check for: if not A, C, G, T, ., *?
    #         if words[4] == 'I' or words[5] == 'I' or words[6] == 'I' or words[7] == 'I' or \
    #             words[4] == 'D' or words[5] == 'D' or words[6] == 'D' or words[7] == 'D':
    #             count += 1
    #             string = words[0] + '\t' + words[2] + '\n'
    #             snp_list.append(string)

    # print('Removing {snps} snps'.format(snps = count))
    # # create text file of positions
    # fileSNP = open('beforePhasing/snpPositionsToRemoveChr' + chr + '.txt', 'w')
    # fileSNP.writelines(snp_list)
    # fileSNP.close()
    # # why dont I exlude-positions? maybe a time saver. currently does 158 seconds
    # print('Creating new target file removed of bad alleles')
    # runCommandLine(['vcftools --vcf beforePhasing/tar/chr' + chr + '.noDups.vcf \
    #     --exclude-positions beforePhasing/snpPositionsToRemoveChr' + chr + '.txt \
    #     --recode \
    #     --out beforePhasing/tar/chr' + chr])