#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 08:47:07 2022
@author: Frederic Labbe
SNP measurement model: This script sample a number of missing loci and randomly replace them with missing data.
It uses a distribution of the proportion of missing SNP loci per host based on Taq-Man assay empirical data.
Note: The distribution file should only contains two columns, i.e. the proportion and the weigths of missing SNP loci per host.
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.txt").
usage: python SnpErrMod.py --inputfile '/path/to/file.txt' --missdist '/path/to/file.txt'
"""

from numpy.random import choice
from random import sample
import os.path
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file in THE REAL McCOIL format.')
parser.add_argument('-m', "--missdist", required = True, help = 'Path to the file providing the distribution of the proportion of missing loci per host.')
args = parser.parse_args()

def SnpErrMod(inputfile, missdist):
    os.chdir(os.path.dirname(os.path.abspath(inputfile)))
    with open(missdist, 'r') as dist:
        header = dist.readline()
        propmiss = list()
        weigth = list()
        for line in dist:
            x = line.split()
            propmiss.append(float(x[0]))
            weigth.append(float(x[1]))
    
    outputfile = inputfile.split("/")[-1].split(".")[0] + "_GenErr" + ".txt"
    f = open(outputfile, 'w')
    with open(inputfile, 'r') as file:
        header = file.readline()
        header = header.split()
        f.write("{}\n".format('\t'.join(map(str, header))))
        header.remove('ind_name')
        nbloci = len(header)
        loci = list(range(1, nbloci + 1))    
        for line in file:
            calls = line.split()
            host = calls[0]
            del calls[0]
            nbmiss = calls.count('-1')
            draw = float(choice(propmiss, 1, p = weigth))
            minmiss = round(draw * nbloci)
            if minmiss > nbmiss:
                difmiss = minmiss - nbmiss
                indexes = sample(loci, difmiss)
                for index in indexes:
                    calls[index -1] = -1
            f.write("{}\t{}\n".format(host, '\t'.join(map(str, calls))))
    f.close()

if __name__ == '__main__':
     SnpErrMod(args.inputfile, args.missdist)