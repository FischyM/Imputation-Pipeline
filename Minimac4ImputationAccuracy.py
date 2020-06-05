#!/usr/bin/enc python3
# Copyright Mathew Fischbach, Myers Lab http://csbio.cs.umn.edu/, 2020
# Code based on the work of Xiaotong Liu from Myers Lab

### Purpose:
# This program will take GWAS vcf files and impute them with Minimac4. An analysis of the
# imputation accuracy is performed.

### Requirements: 
# Input: vcf containing target data and vcf reference panel, both GRCh38 build
# Output: analysis_output.csv (tab-delimited) and per chromosome imputed files in output directory

import sys
import argparse
import os
import subprocess
import json
import csv
from collections import OrderedDict
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


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
    if not os.path.isdir('output'):
        print('directory "output" not found. Creating directory')
        os.mkdir('output')
        for i in range(1, 23):
            os.mkdir('output/chr' + str(i))
        
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
    startDir = os.getcwd()
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
        --prefix afterPhasing/ref/chr' + chr + ' --cpus ' + cpu + ' --rounds 0 --chr ' + chr]
    runCommandLine(toCommandLine)
    print('Done converting ref file to minimac3 format')
    # outputs as .m3vcf.gz

def impute(chr, cpu):
    print('\nBeginning Minimac4 Imputation')
    toCommandLine = ['minimac4 \
        --haps afterPhasing/tar/prePhasedChr' + chr + '.vcf.gz \
        --refHaps afterPhasing/ref/chr' + chr + '.m3vcf.gz \
        --prefix output/chr' + chr + ' \
        --referenceEstimates OFF \
        --mapFile map/genetic_map_hg38_withX.minimac.txt \
        --cpus ' + cpu]
    runCommandLine(toCommandLine)

def setUpCompare(infoFile):
    df = pd.read_csv(infoFile, header=0, sep='\t')
    df[['AvgCall','Rsq', 'LooRsq', 'EmpR', 'EmpRsq', 'Dose0', 'Dose1']] = \
        df[['AvgCall', 'Rsq', 'LooRsq', 'EmpR', 'EmpRsq', 'Dose0', 'Dose1']].apply(pd.to_numeric, \
        errors='coerce')
    return df

def printRsqEmpRsq(df, countGeno, dict_out, maf):
    dict_out['genotyped-snps'][maf]['Rsq'] = round(df['Rsq'].mean(), 4)
    dict_out['genotyped-snps'][maf]['EmpRsq'] = round(df['EmpRsq'].mean(), 4)
    dict_out['genotyped-snps'][maf]['Total snps'] = len(df.index)
    ratio = len(df.index)/float(countGeno)
    dict_out['genotyped-snps'][maf]['Ratio of snps'] = round(ratio, 4)

def rsqCutoffs(mafList, dict_out, type):
    totalSNPs = len(mafList.index)
    cutoffs = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.96, 0.97, 0.98, 0.99]
    for val in cutoffs:
        newMaf = mafList.loc[(mafList['Rsq'] >= val)]
        newSnps = len(newMaf.index)
        dict_out[type]['snp ratios at Rsq cutoff'][round(val,2)] = round(newSnps/float(totalSNPs),4)

def printRsq(df, countTotal, dict_out, maf):
    dict_out['all-snps'][maf]['Rsq'] = round(df['Rsq'].mean(), 4)
    dict_out['all-snps'][maf]['Total snps'] = len(df.index)
    ratio = len(df.index)/float(countTotal)
    dict_out['all-snps'][maf]['Ratio of snps'] = round(ratio, 4)


def compare(infoFile, chr):
    dict_out = {'genotyped-snps': {'MAF <= 0.05%': {}, 'MAF 0.05-5%': {}, \
        'MAF >= 5%': {}, 'correlation': None, 'snp ratios at Rsq cutoff': OrderedDict()}, \
            'all-snps': {'MAF <= 0.05%': {}, 'MAF 0.05-5%': {}, 'MAF >= 5%': {}, \
            'snp ratios at Rsq cutoff': OrderedDict()}}
    df = setUpCompare(infoFile)
    countTotal = len(df.index)
    genotypes = df.loc[(df['Genotyped'] == 'Genotyped')]
    countGeno = len(genotypes.index)

    #1.a
    mafLess_0_05 = df.loc[(df['MAF'] <= 0.0005) & (df['Genotyped'] == 'Genotyped')]
    maf_0_05_to_5 = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] < 0.05) & (df['Genotyped'] == 'Genotyped')]
    mafGreater_5 = df.loc[(df['MAF'] >= 0.05) & (df['Genotyped'] == 'Genotyped')]

    printRsqEmpRsq(mafLess_0_05, countGeno, dict_out, 'MAF <= 0.05%')
    printRsqEmpRsq(maf_0_05_to_5, countGeno, dict_out, 'MAF 0.05-5%')
    printRsqEmpRsq(mafGreater_5, countGeno, dict_out, 'MAF >= 5%')

    #1.b
    corr = mafGreater_5[['Rsq', 'EmpRsq']].corr()
    dict_out['genotyped-snps']['correlation'] = round(corr.iloc[0,1], 4)

    mafGreater_5.plot.scatter(x='Rsq', y='EmpRsq', c='Red')
    plt.savefig('output/chr' + chr + '/scatter_Genotypes.EmpRsq.Rsq_MAF_GT_5.png')

    rsqCutoffs(mafGreater_5, dict_out, 'genotyped-snps')

    #2
    AllMafLess_0_05 = df.loc[(df['MAF'] <= 0.0005)]
    printRsq(AllMafLess_0_05, countTotal, dict_out,'MAF <= 0.05%')

    AllMaf_0_05_to_5 = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] < 0.05)]
    printRsq(AllMaf_0_05_to_5, countTotal, dict_out, 'MAF 0.05-5%')

    AllMafGreater_5 = df.loc[(df['MAF'] >= 0.05)]
    printRsq(AllMafGreater_5, countTotal, dict_out, 'MAF >= 5%')

    #3
    AllMafGreater_5[['Rsq']].plot(kind='hist', bins=100)
    plt.savefig('output/chr' + chr + '/hist_ALL.Rsq_MAF_GT_5.png')

    rsqCutoffs(AllMafGreater_5, dict_out, 'all-snps')

    return dict_out

def writeTSV(dictIn):
    newDict = OrderedDict(dictIn)
    with (open('output/analysis_output.csv', 'w')) as f:
        writer = csv.writer(f, delimiter='\t')
        for a, b in sorted(newDict.items()):
            for i, j in sorted(b.items(), reverse=True):
                for k, l in sorted(j.items()):
                    tempList = []
                    if type(l) is not str and type(l) is not float:
                        for key, val in l.items():
                            tempList.append(str(key) + ':' + str(val))
                        writer.writerow([a, i, k] + tempList)
                    else:
                        writer.writerow([a, i, k, l])


######################################################################################
def main():
    # Set up argument inputs
    targetFileName = ''
    refFileName = ''
    cpus = 1
    analysisDict = {}
    targetFileName, refFileName, cpus = setUp(sys.argv[1:], targetFileName, refFileName, cpus)
    checkPrograms()
    # separateFile(targetFileName, 'tar')
    # separateFile(refFileName, 'ref')
    
    for i in range(1, 23):
        chromosome = str(i)
        # removeDuplicates(chromosome)
        # removeBadAllelesWithAWK(chromosome)
        # prePhase(chromosome, cpus)
        # convertRef(chromosome, cpus)
        # impute(chromosome, cpus)
        infoFile = 'output/chr' + chromosome + '/imputed_chr' + chromosome + '.info'
        if os.path.isfile(infoFile):
            analysisDict['chr' + chromosome] = compare(infoFile, chromosome)
        writeTSV(analysisDict)
if __name__ == '__main__':
    main()

# current testing input
# python Minimac4ImputationAccuracy.py -t PD_NGRC_phs000196.GRCh38.chr1.vcf.gz -r ALL.chr1.shapeit2_integrated_v1a.GRCh38.20181129.phased.ver4.2.vcf.gz --cpus 10
