#!/usr/bin/python
# Copyright Mathew Fischbach, Myers Lab http://csbio.cs.umn.edu/, 2020
# code based heavily on the work of Xiaotong Liu from Myers Lab

# Purpose:
# This program will take plink files and impute them with Minimac4.
# Analysis of the imputation will be done by minimac_info_analysis.py
# This program imputes on a single chromosome at a time. A simple bash script can
# be created to impute chromosomes 1 through 22

# Requirements: 
# Input: vcf file containing target data and vcf reference panel

import sys, getopt
import argparse
import os
import subprocess
import random


startDir = os.getcwd()

def gzipIt(NameOfFile):
    toCommandLine = "gzip -k " + NameOfFile
    runCommandLine([toCommandLine])
    return NameOfFile + '.gz'

def usage():
    print('Minimac4ImputationAccuracy.py -t <target file> -r <ref file> \
        -chr <chromosome> -cpu <threads>\n')

def setUp(argv, targetFileName, refFileName, chromosome, cpus):
    print('Checking input')

    parser = argparse.ArgumentParser(description='Minimac4ImputationAccuracy.py \
        -t <target file> -r <ref file> -chr <chromosome> -cpu <threads>\n')
    parser.add_argument('-t', '--tar', required=True, help='target GWAS dataset as vcf.gz')
    parser.add_argument('-r','--ref', required=True, help='reference dataset as vcf.gz')
    parser.add_argument('-c', '--chr', required=True, help='chromosome to analyze')
    parser.add_argument('--cpus', nargs='?', dest='threads', default=4, help='number of threads to use')
    args = parser.parse_args()
    
    print('Successful input\n')
    targetFileName = args.tar
    refFileName = args.ref
    chromosome = args.chr
    cpus = args.threads
    return targetFileName, refFileName, chromosome, cpus

def runCommandLine(commands):
    try:
        process = subprocess.Popen(commands, shell=True)
        process.communicate()
    except subprocess.CalledProcessError as err:
        print('Error: ', err)
        sys.exit(1)

def checkPrograms():
    print('Checking required programs and directories')
    if not os.path.isdir('vcftools'):
        print('directory "vcftool" not found')
        print('please download or update directory to "vcftools" before proceeding')
        sys.exit(1)

    if not os.path.isdir('Minimac3'):
        print('directory "Minimac3" not found')
        print('please download or update directory to "Minimac3" before proceeding')
        sys.exit(1)
        
    if not os.path.isdir('Minimac4'):
        print('directory "Minimac4" not found')
        print('please download or update directory to "Minimac4" before proceeding')
        sys.exit(1)

    if not os.path.isfile('beagle.27Apr20.b81.jar'):
        print('program  "beagle.27Apr20.b81.jar" not found')
        print('please download "beagle.27Apr20.b81.jar" before proceeding')
        sys.exit(1)

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

    if not os.path.isdir('vcflib'):
        print('directory "vcflib" not found')
        print('please download or update directory to "vcflib" before proceeding')
        sys.exit(1)

    if not os.path.isdir('map'):
        print('directory "map" not found')
        print('please download or update directory to "map" before proceeding. \
            Map should have the needed plink map files which can be downloaded from \
            http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/')
        sys.exit(1)

    print('Successful, all programs/directories are present\n')



def separateFile(file, fileType, chr):
    if file[-2:] != 'gz':
        gz = ''
    else:
        gz = 'gz'
    toCommandLine = ['vcftools/bin/vcftools --' + gz + 'vcf ' + file + ' --chr ' + chr + \
        ' --recode --out beforePhasing/' + fileType + '/chr' + chr ]
    runCommandLine(toCommandLine)
    # target file chromosomes reside in beforePhasing/tar/chr#.recode.vcf
    # reference file chromosomes reside in beforePhasing/ref/chr#.recode.vcf
    # vcftools --recode always outputs as .vcf, no gzipping

def removeDuplicates(chr): # remove duplicate entries. Necessary for pre-phasing and should be in vcf format
    print('Beginning to remove duplicate entries')
    runCommandLine(['vcflib/bin/vcfuniq beforePhasing/tar/chr' + chr + '.recode.vcf > \
        beforePhasing/tar/chr' + chr + '.noDups.vcf'])
    print('Done, duplicate entries removed\n')
    # here I could gzip it for later
    # out: beforePhasing/tar/chr#.noDups.vcf

def removeBadAlleles(chr): 
    # takes snps from target and ref file and creates a new ref file from 1000Genome
    print('finding list of snps from given target and reference files')
    toCommandLine = ['vcftools/bin/vcftools \
        --vcf beforePhasing/tar/chr' + chr + '.noDups.vcf \
        --diff beforePhasing/tar/chr' + chr + '.noDups.vcf \
        --diff-site --out beforePhasing/snpPositionChr' + chr]
    runCommandLine(toCommandLine) 

    # create list of snps while removing poorly coded alleles 
    # this is based off of personal experience. This can easily be changed
    print('\nRemoving wrongly coded alleles')
    snp_list = []
    count = 0
    with open('beforePhasing/snpPositionChr' + chr + '.diff.sites_in_files') as f:
        read_data = f.readlines()
        for i in range(1, len(read_data)): # skip header line
            words = read_data[i].split('\t')
            # Change the logic to check for: if not A, C, G, T, ., *?
            if words[4] == 'I' or words[5] == 'I' or words[6] == 'I' or words[7] == 'I' or \
                words[4] == 'D' or words[5] == 'D' or words[6] == 'D' or words[7] == 'D':
                count += 1
                continue
            if words[1] == '.':
                string = words[0] + '\t' + words[2] + '\n'
                snp_list.append(string)
            else:
                string = words[0] + '\t' + words[1] + '\n'
                snp_list.append(string)

    print('Removing {snps} snps'.format(snps = count))
    # create text file of positions
    fileSNP = open('beforePhasing/snpPositionChr' + chr + '.txt', 'w')
    fileSNP.writelines(snp_list)
    fileSNP.close()
    # why dont I exlude-positions? maybe a time saver. currently does 158 seconds
    print('Creating new target file removed of bad alleles')
    runCommandLine(['vcftools/bin/vcftools --vcf beforePhasing/tar/chr' + chr + '.noDups.vcf \
        --positions beforePhasing/snpPositionChr' + chr + '.txt \
        --recode \
        --out beforePhasing/tar/chr_ready_to_prephase' + chr])

def prePhase(chr, cpu):
    print('beginning prephasing with beagle')
    toCommandLine = ['java -jar beagle.27Apr20.b81.jar gt=beforePhasing/tar/chr_ready_to_prephase' + chr + '.recode.vcf \
        out=afterPhasing/tar/prePhasedChr' + chr + \
        ' iterations=20 map=map/plink.chr' + chr + '.GRCh38.map \
        chrom=' + chr + \
        ' nthreads=' + cpu]
    runCommandLine(toCommandLine)
    print('Done prephasing')
    # outputs as .vcf.gz

def convertRef(chr, cpu):
    # to run minimac4, ref file must be in m3vcf format. uses minimac3 to do this.
    # It was determined that we will let Minimac3 do its parameter estimation so that a map
    # file input is no longer needed. This makes the process easier when dealing with
    # different genome builds
    print('\nBeginning reference file conversion to m3vcf')
    toCommandLine = ['Minimac3/bin/Minimac3 --refHaps beforePhasing/ref/chr' + chr + '.recode.vcf \
        --processReference \
        --prefix afterPhasing/ref/chr' + chr + ' --cpus ' + cpu + ' --rounds 5 --chr ' + chr]
    runCommandLine(toCommandLine)
    print('Done converting ref file to minimac3 format')
    # outputs as .m3vcf.gz

def impute(chr, cpu):
    # shouldn't need --ignore-duplicates since that was taken care of
    toCommandLine = ['Minimac4/release-build/minimac4 --haps afterPhasing/tar/prePhasedChr' \
        + chr + '.vcf.gz \
        --refHaps afterPhasing/ref/chr' + chr + '.m3vcf.gz \
        --prefix afterPhasing/imputed_chr' + chr + \
        ' --cpus ' + cpu + ' --chr ' + chr]
    runCommandLine(toCommandLine)
    print('imputation using minimac4 - done')

def compare(chr):
    os.chdir('afterPhasing')
    runCommandLine(['python3 ../minimac_info_analysis.py afterPhasing/imputed_chr' + chr + '.info \
        > imputation_assessment_chr' + chr + '.txt'])
    # simply sends stdout to text file 
    os.chdir(startDir)



######################################################################################
if __name__ == '__main__':
    # Set up argument inputs
    targetFileName = ''
    refFileName = ''
    chromosome = ''
    cpus = ''
    targetFileName, refFileName, chromosome, cpus = \
        setUp(sys.argv[1:], targetFileName, refFileName, chromosome, cpus)

    # Check for required programs and directories 
    checkPrograms()

    #separateFile(targetFileName, 'tar', chromosome)
    #separateFile(refFileName, 'ref', chromosome)

    removeDuplicates(chromosome)
    removeBadAlleles(chromosome)
    prePhase(chromosome, cpus)
    # convertRef(chromosome, cpus)
    # impute(chromosome, cpus)
    # compare(chromosome)






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