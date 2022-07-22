#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 09:07:12 2022
@author: Frederic Labbe
This script calculates the proportions of calls per host and per locus from a file in THE REAL McCOIL format.
Two output files will be generated, i.e.:
-"file_name_PropCalls_Ind.txt": proportions of calls per host,
-"file_name_PropCalls_Loc.txt": proportions of calls per locus.
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.txt").
usage: python CalculPropCalls.py --inputfile '/path/to/file.txt'
"""

import os.path
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file in THE REAL McCOIL format.')
args = parser.parse_args()

def PropCalls(inputfile):
    if os.path.exists(inputfile):
        os.chdir(os.path.dirname(os.path.abspath(inputfile)))
        outputfile1 = inputfile.split("/")[-1].split(".")[0] + "_PropCalls_Ind" + ".txt"
        f1 = open(outputfile1, 'w')
        with open(inputfile, 'r') as file:
            header = file.readline()
            header = header.split()
            header.remove('ind_name')
            nbloci = len(header)
            listloci = [[] for _ in range(nbloci)]
            f1.write("ind_name\tprop_miss\tprop_homo_min\tprop_hetero\tprop_homo_maj\n")
            for line in file:
                calls = line.split()
                host = calls[0]
                del calls[0]
                miss = calls.count('-1') / len(calls)
                homo_min = calls.count('0') / len(calls)
                homo_maj = calls.count('1') / len(calls)
                het = calls.count('0.5') / len(calls)
                f1.write("{}\t{}\t{}\t{}\t{}\n".format(host, miss, homo_min, het, homo_maj))
                for locus in range(0, nbloci, 1):
                    listloci[locus].append(calls[locus])
        f1.close()
        outputfile2 = inputfile.split("/")[-1].split(".")[0] + "_PropCalls_Loc" + ".txt"
        f2 = open(outputfile2, 'w')
        f2.write("loc_name\tprop_miss\tprop_homo_min\tprop_hetero\tprop_homo_maj\n")
        for locus in range(0, nbloci, 1):
                miss = listloci[locus].count('-1') / len(listloci[locus])
                homo_min = listloci[locus].count('0') / len(listloci[locus])
                homo_maj = listloci[locus].count('1') / len(listloci[locus])
                het = listloci[locus].count('0.5') / len(listloci[locus])
                f2.write("{}\t{}\t{}\t{}\t{}\n".format(locus + 1, miss, homo_min, het, homo_maj))
        f2.close()
    else:
       print('Error: provide a valid path to the input file')

if __name__ == '__main__':
     PropCalls(args.inputfile)