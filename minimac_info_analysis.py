
# Copyright Mathew Fischbach, Myers Lab http://csbio.cs.umn.edu/, 2020
# Purpose: analyze the .info file that results from a Minimac4 imputation
# Requirements: pandas, numpy, matplotlib, python3
import json
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def setUp(infoFile):
    # infoFile = sys.argv[1] # take in the .info file
    df = pd.read_csv(infoFile, header=0, sep='\t')

    df[['AvgCall','Rsq', 'LooRsq', 'EmpR', 'EmpRsq', 'Dose0', 'Dose1']] = \
        df[['AvgCall', 'Rsq', 'LooRsq', 'EmpR', 'EmpRsq', 'Dose0', 'Dose1']].apply(pd.to_numeric, \
        errors='coerce')

    return df

def printRsqEmpRsq(df, countGeno, dict_out, maf):
    dict_out['genotyped snps'][maf]['Rsq'] = round(df['Rsq'].mean(), 4)
    dict_out['genotyped snps'][maf]['EmpRsq'] = format(round(df['EmpRsq'].mean(), 4))
    dict_out['genotyped snps'][maf]['Total snps'] = format(len(df.index))
    ratio = float(len(df.index)/countGeno)
    dict_out['genotyped snps'][maf]['Ratio of snps'] = format(round(ratio, 4))

def rsqCutoffs(mafList, dict_out, type):
    totalSNps = len(mafList.index)
    cutoffs = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.96, 0.97, 0.98, 0.99]
    for val in cutoffs:
        newMaf = mafList.loc[(mafList['Rsq'] >= val)]
        newSnps = len(newMaf.index)
        dict_out[type]['snp ratios at Rsq cutoff'][round(val,2)] = round((newSnps/totalSNps),4)

def printRsq(df, coutnTotal, dict_out, maf):
    dict_out['all snps'][maf]['Rsq'] = round(df['Rsq'].mean(), 4)
    dict_out['all snps'][maf]['Total snps'] = len(df.index)
    ratio = float(len(df.index)/coutnTotal)
    dict_out['all snps'][maf]['Ratio of snps'] = round(ratio, 4)


def run_compare(infoFile):
    dict_out = {'genotyped snps': {'MAF <= 0.05%': {}, 'MAF > 0.05% and < 5%': {}, \
        'MAF >= 5%': {}, 'correlation': None, 'snp ratios at Rsq cutoff': {}}, \
            'all snps': {'MAF <= 0.05%': {}, 'MAF > 0.05% and < 5%': {}, 'MAF >= 5%': {}, \
            'snp ratios at Rsq cutoff': {}}}
    df = setUp(infoFile)
    coutnTotal = len(df.index)
    genotypes = df.loc[(df['Genotyped'] == 'Genotyped')]
    countGeno = len(genotypes.index)

    #1.a
    mafLess_0_05 = df.loc[(df['MAF'] <= 0.0005) & (df['Genotyped'] == 'Genotyped')]
    maf_0_05_to_5 = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] < 0.05) & (df['Genotyped'] == 'Genotyped')]
    mafGreater_5 = df.loc[(df['MAF'] >= 0.05) & (df['Genotyped'] == 'Genotyped')]

    printRsqEmpRsq(mafLess_0_05, countGeno, dict_out, 'MAF <= 0.05%')
    printRsqEmpRsq(maf_0_05_to_5, countGeno, dict_out, 'MAF > 0.05% and < 5%')
    printRsqEmpRsq(mafGreater_5, countGeno, dict_out, 'MAF >= 5%')

    #1.b
    corr = mafGreater_5[['Rsq', 'EmpRsq']].corr()
    dict_out['genotyped snps']['correlation'] = round(corr.iloc[0,1], 4)

    mafGreater_5.plot.scatter(x='Rsq', y='EmpRsq', c='Red')
    plt.savefig('scatter_Genotypes.EmpRsq.Rsq_MAF_GT_5.png')

    rsqCutoffs(mafGreater_5, dict_out, 'genotyped snps')

    #2
    AllMafLess_0_05 = df.loc[(df['MAF'] <= 0.0005)]
    printRsq(AllMafLess_0_05, coutnTotal, dict_out,'MAF <= 0.05%')

    AllMaf_0_05_to_5 = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] < 0.05)]
    printRsq(AllMaf_0_05_to_5, coutnTotal, dict_out, 'MAF > 0.05% and < 5%')

    AllMafGreater_5 = df.loc[(df['MAF'] >= 0.05)]
    printRsq(AllMafGreater_5, coutnTotal, dict_out, 'MAF >= 5%')

    #3
    AllMafGreater_5[['Rsq']].plot(kind='hist', bins=100)
    plt.savefig('hist_ALL.Rsq_MAF_GT_5.png')

    rsqCutoffs(AllMafGreater_5, dict_out, 'all snps')

    return dict_out
