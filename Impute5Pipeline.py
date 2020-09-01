#!/usr/bin/env python3.8
# # Copyright Mathew Fischbach, Myers Lab http://csbio.cs.umn.edu/, 2020
# Code based on the work of Xiaotong Liu from Myers Lab

# Purpose:
# This program imputes GWAS data files.
# An analysis of the imputation accuracy is also performed.

# Supporting Information:
# 1000 Genome GRCh38 reference files with RSIDs: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/

from functools import partial
from multiprocessing.dummy import Pool

import allel
import sys
import argparse
import os
import subprocess
import csv
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt


def setUp(argv, targetFileName, refFileName, cpus, prefix, mem):
    print('\nChecking input')
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-t', '--tarHaps', type=str, required=True, help='Target GWAS dataset in vcfn vcf bgziped, or bcf format as a singular file. This program will do the separation')
    parser.add_argument(
        '-r', '--refHaps', type=str, required=True, help='File list containing reference datasets (in order of chromosomes 1-22, one per line) in imp5 format.')
    parser.add_argument(
        '-p', '--prefix', type=str, default='ouput', help='Prefix for output files.')
    parser.add_argument(
        '-c', '--cpus', type=str, default='1', help='Number of threads to use. Suggested around 12 cores with approx 100 Gb RAM.')
    parser.add_argument(
        '-m', '--memory', type=str, default='10000000', help='length of region to impute at a time. 10 cM is default. To reduce memory usage, use less than 10000000.')
    args = parser.parse_args()

    targetFileName = os.path.abspath(args.tarHaps)
    refFileName = os.path.abspath(args.refHaps)
    cpus = args.cpus
    if int(cpus) < 4:
        print('Warning: the number of cpus to use is less than 4. \
        This will be very limiting to this pipeline.\nProgram will continue \
        to run unless explicitly stopped.')
    cpus = args.cpus
    prefix = args.prefix
    try:
        int(args.memory)
    except:
        print('please supply a number to -m/--memory')
        sys.exit(1)
    mem = args.memory
    print('Successful\n')

    return targetFileName, refFileName, cpus, prefix, mem


def runCommandLine(commands):
    try:
        subprocess.run(commands, shell=True, check=True)
    except subprocess.CalledProcessError as err:
        print('Error: ', err)
        sys.exit(1)

def checkPrograms():
    print('Checking required programs and directories')

    # create directories if not present
    if not os.path.isdir('output'):
        print('directory "output" not found. Creating temporary directory')
        os.mkdir('output')
        for i in range(1, 23):
            os.mkdir('output/chr' + str(i))

    if not os.path.isdir('beforePhasing'):
        print('directory "beforePhasing" not found. Creating temporary directory')
        os.mkdir('beforePhasing')

    if not os.path.isdir('afterPhasing'):
        print('directory "afterPhasing" not found. Creating temporary directory')
        os.mkdir('afterPhasing')

    if not os.path.isdir('map'):
        print('directory for map data not found, be sure to download from the git repository, unzip plink.GRCh23.map.zip and have the files reside in map/')
        sys.exit(1)

    # set all programs used to the env PATH variable
    if 'htslib' not in os.environ.get('PATH'):
        print('Cannot find htslib on PATH, please add it. Bgzip and tabix are used with this.')
        sys.exit(1)

    if 'bcftools' not in os.environ.get('PATH'):
        print('Cannot find bcftools on PATH, please add it. Used for sorting and removing duplicates')
        sys.exit(1)
    
    if 'impute_v5.1' not in os.environ.get('PATH'):
        print('Cannot find impute5 on PATH, please add it. Used for imputation')
        sys.exit(1)

    if 'shapeit4' not in os.environ.get('PATH'):
        print('Cannot find shapeit4 on PATH, please add it. Used for prephasing')
        sys.exit(1)

    print('Successful\n')


def splitVCF(file, cpus):
    # as long as the target fiile is vcf, vcf.gz. or bcf, bcftools sort will run and convert to bcf
    if not os.path.isfile('beforePhasing/change_chromosomes_tmp.txt'):
        with open('beforePhasing/change_chromosomes_tmp.txt', 'w') as f:
            for i in range(1, 23):
                line = str(i) + " chr" + str(i) + "\n"
                f.write(line)
    
    runCommandLine(["bcftools sort -Ou -T beforePhasing/ " + file + 
    " | bcftools annotate --rename-chrs beforePhasing/change_chromosomes_tmp.txt -Ob -o beforePhasing/sortedTargetFile.tmp.bcf --threads " + cpus])
    runCommandLine(['tabix -f beforePhasing/sortedTargetFile.tmp.bcf'])

    print('separating chromosomes 1-22')
    for i in range(1, 23):
        chrm = str(i)
        runCommandLine(["bcftools view --threads " + cpus + " -Ob -o beforePhasing/chr" + chrm + ".tmp.bcf \
            -r chr" + chrm + " beforePhasing/sortedTargetFile.tmp.bcf"])

def sortAndFilterShapeit4(chrm, cpus):
    # remove duplicate. File should be vcf.gz
    print('Checking QC: removing duplicates, aligning REF and ALT alleles')
    runCommandLine(['bcftools norm beforePhasing/chr' + chrm + '.tmp.bcf -d snps -Ob -c s -N \
        -f /project/csbio/MathewF/imputation/ImputationAccuracy/reference_1000G/1000G_fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa \
        --threads ' + cpus + " -o beforePhasing/chr" + chrm + ".prephaseReady.tmp.bcf"])
    runCommandLine(["tabix beforePhasing/chr" + chrm + ".prephaseReady.tmp.bcf"])
    # fasta alignment file here. Added above to bcftools norm command
    # http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

def prePhaseShapeit4(chrm, cpu):
    print('beginning prephasing with Shapeit4')
    runCommandLine(["shapeit4 --input beforePhasing/chr" + chrm + ".prephaseReady.tmp.bcf \
        --map map/shapeit4_maps/chr" + chrm + ".b38.gmap.gz \
        --region chr" + chrm + " \
        --thread " + cpu + " \
        --output afterPhasing/prePhasedChr" + chrm + "_needGENO.tmp.bcf --log afterPhasing/prePhasedChr" + chrm + ".tmp.log"])

def runImpute5(chrm, refFile, cpus, mem):
    print('\nBeginning Impute5 Imputation')
    region = int(mem)  # 10 cM default

    # Find how long the chromosome is using the GRCh38 fasta file index
    commandHeader = "grep -w chr" + chrm + \
    " reference_1000G/1000G_fasta/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai | cut -f 2"
    try:
        process = subprocess.Popen(commandHeader, shell=True, stdout=subprocess.PIPE, text=True, universal_newlines=True)
    except subprocess.CalledProcessError as err:
        print('Error: ', err)
        sys.exit(1)
    length = process.communicate()
    length = int(length[0])
    div = abs(-(length) // region)

    # need index for target file
    runCommandLine(["tabix -f afterPhasing/prePhasedChr" + chrm + ".tmp.bcf"])

    commands = []
    imputeList = []
    for i in range(0, div):
        start = 1 if i == 0 else (i * region) + 1
        end = length if i == (int(div) - 1) else (i + 1) * region
        step = str(i + 1)
        commandLine = ["impute5 --g afterPhasing/prePhasedChr" + chrm + ".tmp.bcf \
            --h " + refFile.strip() + " --m map/impute5_maps/chr" + chrm + ".b38.gmap.gz \
            --r chr" + chrm + ":" + str(start) + "-" + str(end) + " --o output/chr" + 
            chrm + "/impute5_" + step + ".tmp.bcf \
            --l output/chr" + chrm + "/impute5_" + step + ".log --out-gp-field"]
        commands.append(commandLine)
        imputeList.append("output/chr" + chrm + "/impute5_" + step + ".tmp.bcf")
    # the following code was found at https://stackoverflow.com/questions/14533458/python-threading-multiple-bash-subprocesses/14533902#14533902
    max_workers = int(cpus)
    pool = Pool(max_workers) # run the value of 'cpus' concurrent commands at a time
    for i, returncode in enumerate(pool.imap(partial(subprocess.call, shell=True, stdout=subprocess.DEVNULL), commands)):
        if returncode != 0:
            print("region %d failed: %d\nView output/Chr%d/impute5_%d.log for more details" % ((i+1), returncode, int(chrm), (i+1)))

    # concat vcf files together
    imputeString = ''
    for file in imputeList:
        if os.path.isfile(file):
            imputeString = imputeString + file + " "
    runCommandLine(["bcftools concat --no-version -Ob --threads " + cpus + 
        " -o output/chr" + chrm + "/impute5_all.tmp.bcf " + imputeString])


def setUpCompareImpute5(cpus, chrm):
    infoFile = "output/chr" + chrm + "/impute5_all_MAF_GENO.tmp.vcf.gz"
    # bcftools MAF calculation 25:69min. Can bring it down to 11min if I filter by MAF without using fill-tags
    runCommandLine(["bcftools +fill-tags --no-version --threads " + cpus + " output/chr" + chrm + "/impute5_all.tmp.bcf -Ob -o output/chr" + chrm + "/impute5_all_MAF.tmp.bcf -- -t MAF"])
    
    # find which variants are genotyped 7:7min
    runCommandLine(["bcftools query -f '%CHROM\t%POS\t1\n' -o output/chr" + chrm + "/region_list_tmp.txt afterPhasing/prePhasedChr" + chrm + "_needGENO.tmp.bcf"])
    runCommandLine(["bgzip -f output/chr" + chrm + "/region_list_tmp.txt; tabix -s1 -b2 -e2 output/chr" + chrm + "/region_list_tmp.txt.gz"])
    # Create annotation file for markers that are similar b/t before and after imputation, and filter for IMP flag
    runCommandLine(["bcftools query -T output/chr" + chrm + "/region_list_tmp.txt.gz output/chr" + chrm + "/impute5_all_MAF.tmp.bcf -f '%CHROM\t%POS\t1\n' -o output/chr" + chrm + "/annot_geno_list.tmp.txt"])
    runCommandLine(["bgzip -f output/chr" + chrm + "/annot_geno_list.tmp.txt; tabix -s1 -b2 -e2 output/chr" + chrm + "/annot_geno_list.tmp.txt.gz"])

    # Change header lines
    runCommandLine(["""echo '##INFO=<ID=GENO,Number=0,Type=Flag,Description="Genotyped marker">' > output/chr""" + chrm + "/header_lines_tmp"])
    # Annotate imputed file 7:52min
    runCommandLine(["bcftools annotate --no-version -a output/chr" + chrm + "/annot_geno_list.tmp.txt.gz -h output/chr" + chrm + "/header_lines_tmp -c CHROM,POS,INFO/GENO output/chr" + chrm + "/impute5_all_MAF.tmp.bcf \
        -Ob -o output/chr" + chrm + "/impute5_all_MAF_GENO.tmp.bcf --threads " + cpus])
    # prepare INFO like file (just remove the genotypes)1:10min
    runCommandLine(["bcftools view --no-version -G output/chr" + chrm + "/impute5_all_MAF_GENO.tmp.bcf -Oz -o " + infoFile + " --threads " + cpus])
    
    print("\topening info file as a data frame...")
    df = allel.vcf_to_dataframe(infoFile, fields=['variants/*'])
    return df, infoFile

def printInfoImpute(dfInput, typeVar, csvList, maf, chrm):
    csvList.append([chrm, typeVar , maf, round(dfInput['INFO'].mean(), 4)])

def createGraphs(df, infoFile, chrm):
    print('new plot')
    binNum = 30
    all_maf_0_1 = df.loc[(df['MAF'] <= 0.01)]
    all_maf_GT_1 = df.loc[(df['MAF'] > 0.01)]
    bin_means_info_0_1 = stats.binned_statistic(all_maf_0_1['MAF'], all_maf_0_1['INFO'], 'mean', bins=binNum)
    bin_means_info_GT_1 = stats.binned_statistic(all_maf_GT_1['MAF'], all_maf_GT_1['INFO'], 'mean', bins=binNum)
    xaxis_values_0_1  = [ x for x in np.linspace(0.0, 0.01, num=binNum, endpoint=True)]
    xaxis_values_GT_1  = [ x for x in np.linspace(0.01, 0.50, num=binNum, endpoint=True)]


    plt.subplot(121)
    plt.plot(xaxis_values_0_1, bin_means_info_0_1[0], 'r--', label='INFO score')
    plt.title('Average INFO, MAF 0-1%')
    plt.xlabel('log(MAF)')
    plt.ylabel('Average INFO')
    plt.grid(True)
    plt.legend()
    plt.xscale('log')
    plt.tight_layout(pad=1.0)
    plt.subplot(122)
    plt.plot(xaxis_values_GT_1, bin_means_info_GT_1[0], 'r--', label='INFO score')
    plt.title('Average INFO, MAF 1-50%')
    plt.xlabel('log(MAF)')
    plt.ylabel('Average INFO')
    plt.grid(True)
    plt.legend()
    plt.xscale('log')
    plt.tight_layout(pad=1.0)
    plt.savefig(infoFile[:-7] + '_infoScoreByMaf.jpeg')

def compareImpute5(cpus, chrm):
    print("Starting accuracy analysis")
    df, infoFile = setUpCompareImpute5(cpus, chrm)
    df.rename(columns={"MAF_1": "MAF"}, inplace=True)
    print("\tbeginning accuracy calculations...")
    csvList = []

    mafLess_0_05_GENO = df.loc[(df['MAF'] <= 0.0005) & (df['IMP'] == False) & (df['GENO'] == True) ]
    maf_0_05_to_5_GENO = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] <= 0.05) & (df['IMP'] == False) & (df['GENO'] == True) ]
    maf_5_to_50_GENO = df.loc[(df['MAF'] > 0.05) & (df['MAF'] <= 0.5) & (df['IMP'] == False) & (df['GENO'] == True) ]
    printInfoImpute(mafLess_0_05_GENO, 'genotyped', csvList, '0-0.05%', chrm)
    printInfoImpute(maf_0_05_to_5_GENO, 'genotyped', csvList, '0.05-5%', chrm)
    printInfoImpute(maf_5_to_50_GENO, 'genotyped', csvList, '5-50%', chrm)

    mafLess_0_05_IMP = df.loc[(df['MAF'] <= 0.0005) & (df['IMP'] == True) & (df['GENO'] == False) ]
    maf_0_05_to_5_IMP = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] <= 0.05) & (df['IMP'] == True) & (df['GENO'] == False) ]
    maf_5_to_50_IMP = df.loc[(df['MAF'] > 0.05) & (df['MAF'] <= 0.5) & (df['IMP'] == True) & (df['GENO'] == False) ]
    printInfoImpute(mafLess_0_05_IMP, 'imputed', csvList, '0-0.05%', chrm)
    printInfoImpute(maf_0_05_to_5_IMP, 'imputed', csvList, '0.05-5%', chrm)
    printInfoImpute(maf_5_to_50_IMP, 'imputed', csvList, '5-50%', chrm)

    mafLess_0_05_all = df.loc[(df['MAF'] <= 0.0005) ]
    maf_0_05_to_5_all = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] <= 0.05) ]
    maf_5_to_50_all = df.loc[(df['MAF'] > 0.05) & (df['MAF'] <= 0.5) ]
    printInfoImpute(mafLess_0_05_all, 'all', csvList, '0-0.05%', chrm)
    printInfoImpute(maf_0_05_to_5_all, 'all', csvList, '0.05-5%', chrm)
    printInfoImpute(maf_5_to_50_all, 'all', csvList, '5-50%', chrm)

    print("Number of IMP flagged markers {}".format(len(df.loc[(df['IMP'] == True) & (df['GENO'] == False)])))
    print("Number of GENO flagged markers {}".format(len(df.loc[(df['IMP'] == False) & (df['GENO'] == True)])))
    print("Number of IMP and GENO flagged markers {}".format(len(df.loc[(df['IMP'] == True) & (df['GENO'] == True)])))
    print("Number of no IMP and no GENO flagged markers {}".format(len(df.loc[(df['IMP'] == False) & (df['GENO'] == False)])))

    createGraphs(df, infoFile, chrm)

    print("Successful")
    return csvList

def writeCSV_Impute5(analysisList, chrm):
    with (open('output/chr' + chrm + '/Impute5_analysis_output_tmp.csv', 'w')) as f:
        writer = csv.writer(f)
        writer.writerows(analysisList)

def combineCSVsImpute5(prefix, outputDir):
    header = "chr, SNP type, MAF, AvgINFO"
    runCommandLine(['echo ' + header + ' > ' + outputDir + '/' + prefix + '_imputation_accuracy.csv'])
    runCommandLine(['cat output/chr*/Impute5_analysis_output_tmp.csv >> ' + outputDir + '/' + prefix + '_imputation_accuracy.csv'])

def cleanUp(prefix, cpus, outputDir):
    print("concatenating all chromosomes...")
    chrList = ''
    for i in range(1, 23):
        chrList = chrList + "output/chr" + str(i) + "/impute5_all_MAF_GENO.tmp.bcf "
    runCommandLine(["bcftools concat --no-version -Ob --threads " + cpus + " -o " + outputDir + '/' + prefix + ".bcf " + chrList])
    print("removing temporary files and directories...")
        

def main():
    # Set up argument inputs
    targetFileName = ''
    prefix = ''
    refFileName = ''
    cpus = ''
    mem =''
    
    targetFileName, refFileName, cpus, prefix, mem = setUp(sys.argv[1:], targetFileName, refFileName, cpus, prefix, mem)
    
    with open(refFileName) as f:
        refFileList = f.readlines()
    tempDir = os.path.dirname(sys.argv[0])
    outputDir = os.getcwd()
    os.chdir(tempDir)

    checkPrograms()
    splitVCF(targetFileName, cpus)

    # Begin working on each chromosome
    for i in range(1, 23):
        chromosome = str(i)
        print('Starting on chromosome ' + chromosome)

        sortAndFilterShapeit4(chromosome, cpus)
        prePhaseShapeit4(chromosome, cpus)

        refFile = os.path.dirname(refFileName) + "/" + refFileList[i-1]
        runImpute5(chromosome, refFile, cpus, mem)
        print("Successful chromosome " + chromosome + " imputation\n")

        analysisList = compareImpute5(cpus, chromosome)
        writeCSV_Impute5(analysisList, chromosome)
        print('Successful accuracy analysis\n')
    
    combineCSVsImpute5(prefix, outputDir)
    cleanUp(prefix, cpus, outputDir)
    os.chdir(outputDir)

if __name__ == '__main__':
    main()


# ../../ImputationAccuracy/Impute5Pipeline.py -t PD_CIDR_phs000126.GRCh38.vcf.gz -r ../../ImputationAccuracy/reference_1000G/impute5Ref/refFileList -c 12 -p outputTesting

