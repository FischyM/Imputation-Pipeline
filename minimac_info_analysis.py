#!/usr/bin/python
# Copyright Mathew Fischbach, Myers Lab http://csbio.cs.umn.edu/, 2020
# Purpose: analyze the .info file that results from a Minimac4 imputation
# Requirements: pandas, numpy, matplotlib, python3


import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def setUp():
    infoFile = sys.argv[1] # take in the .info file
    df = pd.read_csv(infoFile, header=0, sep='\t')

    print('Retyping column data to float type')
    df[['AvgCall','Rsq', 'LooRsq', 'EmpR', 'EmpRsq', 'Dose0', 'Dose1']] = \
        df[['AvgCall', 'Rsq', 'LooRsq', 'EmpR', 'EmpRsq', 'Dose0', 'Dose1']].apply(pd.to_numeric, \
        errors='coerce')

    return df

def printRsqEmpRsq(df, countGeno):
    print('\tRsq: {}'.format(round(df['Rsq'].mean(), 4)))
    print('\tEmpRsq: {}'.format(round(df['EmpRsq'].mean(), 4)))
    print('\tTotal snps: {}'.format(len(df.index)))
    ratio = float(len(df.index)/countGeno)
    print('\tRatio of snps: {}'.format(round(ratio, 4)))

def rsqCutoffs(mafList):
    totalSNps = len(mafList.index)
    cutoffs = [0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.96, 0.97, 0.98, 0.99]
    for val in cutoffs:
        newMaf = mafList.loc[(mafList['Rsq'] >= val)]
        newSnps = len(newMaf.index)
        print('\tRatio of snps at rsq cutoff of {value}: {total}'.format(value = round(val,2), total = round((newSnps/totalSNps),4)))

def printRsq(df, coutnTotal):
    print('\tRsq: {}'.format(round(df['Rsq'].mean(), 4)))
    print('\tTotal snps: {}'.format(len(df.index)))
    ratio = float(len(df.index)/coutnTotal)
    print('\tRatio of snps: {}'.format(round(ratio, 4)))


def main():
    df = setUp()

    coutnTotal = len(df.index)

    genotypes = df.loc[(df['Genotyped'] == 'Genotyped')]
    countGeno = len(genotypes.index)


    #1.a
    mafLess_0_05 = df.loc[(df['MAF'] <= 0.0005) & (df['Genotyped'] == 'Genotyped')]
    maf_0_05_to_5 = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] < 0.05) & (df['Genotyped'] == 'Genotyped')]
    mafGreater_5 = df.loc[(df['MAF'] >= 0.05) & (df['Genotyped'] == 'Genotyped')]

    print('Looking at Rsq and EmpRsq for genotyped snps')

    print('\nMAF less than or equal to 0.05%')
    printRsqEmpRsq(mafLess_0_05, countGeno)

    print('\nMAF greater than 0.05%, less than 5%')
    printRsqEmpRsq(maf_0_05_to_5, countGeno)

    print('\nMAF greater than or equal to 5%')
    printRsqEmpRsq(mafGreater_5, countGeno)
    print()


    #1.b
    print('Correlation between Rsq and EmpRsq')
    print(mafGreater_5[['Rsq', 'EmpRsq']].corr())
    print()

    mafGreater_5.plot.scatter(x='Rsq', y='EmpRsq', c='Red')
    plt.savefig('scatter_Genotypes.EmpRsq.Rsq_MAF_GT_5.png')

    rsqCutoffs(mafGreater_5)

    #2
    print('\nLooking at Rsq for ALL snps')

    print('\nMAF less than or equal to 0.05%')
    AllMafLess_0_05 = df.loc[(df['MAF'] <= 0.0005)]
    printRsq(AllMafLess_0_05, coutnTotal)

    print('\nMAF greater than 0.05%, less than 5%')
    AllMaf_0_05_to_5 = df.loc[(df['MAF'] > 0.0005) & (df['MAF'] < 0.05)]
    printRsq(AllMaf_0_05_to_5, coutnTotal)

    print('\nMAF greater than or equal to 5%')
    AllMafGreater_5 = df.loc[(df['MAF'] >= 0.05)]
    printRsq(AllMafGreater_5, coutnTotal)
    print()

    #3

    AllMafGreater_5[['Rsq']].plot(kind='hist', bins=100)
    plt.savefig('hist_ALL.Rsq_MAF_GT_5.png')

    rsqCutoffs(AllMafGreater_5)


if __name__ == '__main__':
    main()