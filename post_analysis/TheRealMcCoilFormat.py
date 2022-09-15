#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 10:13:17 2022
@author: flabbe
Based on the SNP allele frequencies, this script converts the SNP data into THE REAL McCOIL input format.
See the README.docx file of THEREALMcCOIL github for more details; link: https://github.com/EPPIcenter/THEREALMcCOIL
SNP calling information is stored in a matrix, where each element Sij represents SNP information at locus j of individual i,
and can be 0 [homozygous minor allele], 0.5 [heterozygous], 1 [homozygous major allele] or -1 [missing data].
It uses the "sampled_infections" and "sampled_infection_snps" tables from the varmodel3 output database (in SQLite3 format).
Note: The script also exports the SNP minor allele frequencies (MAF) in another output file.
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
usage: python TheRealMcCoilFormat.py --inputfile '/path/to/file.txt' --time 300 --minfreq 0.1
"""

import os.path
import sqlite3
import pandas as pd
import argparse
import sys
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file')
parser.add_argument('-t', "--time", type = int, required = True, help = 'Time to make the calculations')
parser.add_argument('-m', "--minfreq", type = float, required = True, help = 'Minor allele frequency (MAF) for a SNP to be considered')
args = parser.parse_args()

def TheRealMcCoilFormat(inputfile, time, minfreq):
    if os.path.exists(inputfile):
        os.chdir(os.path.dirname(os.path.abspath(inputfile)))
        if 0 <= minfreq <= 0.5:

            # Extract and reformat information from the output database
            con = sqlite3.connect(inputfile)
            df1 = pd.read_sql_query('SELECT time, host_id, infection_id, strain_id, expression_index FROM sampled_infections', con)
            df1.columns = ['time', 'host_id', 'infection_id', 'strain_id', 'expression_index_infection']
            df1 = df1.dropna()
            df2 = pd.read_sql_query('SELECT infection_id, snp_index, allele_snp_id FROM sampled_infection_snps', con)
                   
            # Merge and filter
            outputfile1 = inputfile.split("/")[-1].split(".")[0] + "_SNPfreq_" + str(time) + "days.txt"
            f = open(outputfile1, 'w')
            f.write("snp_index\tnb_allele_0\tnb_allele_1\tnb_host\tfreq_allele_0\tfreq_allele_1\tmaf\n")
            df = df1.merge(df2, left_on = 'infection_id', right_on = 'infection_id')
            maf_dict = {}
            if df['time'].isin([time]).any():
                df_time = df[df['time'] == time]
                
                # Calculate the allele frequencies
                header = 'ind_name'
                for index in list(df_time['snp_index'].unique()):
                    snp = df_time[df_time['snp_index'] == index]
                    all1 = list(snp.allele_snp_id).count(1)
                    all2 = list(snp.allele_snp_id).count(2)
                    host = len(list(snp.allele_snp_id))
                    freqall1 = all1 / host
                    freqall2 = 1 - freqall1
                    maf = min(freqall1, freqall2)
                    if  all1 > all2:
                        major = 1
                    else:
                        major = 2
                    maf_dict[index] = maf, major
                    header = header + '\tsite' + str(index)
                    
                    # Export the allele frequencies
                    f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(index, all1, all2, host, freqall1, freqall2, maf))
                    
                # Converte data into THE REAL McCOIL input format
                header = header + '\n'
                outputfile2 = inputfile.split("/")[-1].split(".")[0] + "_TheRealMcCOIL_Input_" + str(time) + "days.txt"
                g = open(outputfile2, 'w')
                g.write(header)
                for ind in sorted(list(df_time['host_id'].unique()), key = int):
                    df_time_ind = df_time[df_time['host_id'] == ind]
                    individual = 'ind' + str(ind)
                    for index in list(df_time['snp_index'].unique()):
                        df_time_ind_index = df_time_ind[df_time_ind['snp_index'] == index]
                        if maf_dict[index][0] < minfreq:
                            individual = individual + '\t-1'
                        else:
                            if len(df_time_ind_index.index) < 2:
                                if maf_dict[index][1] == df_time_ind_index.iloc[0, 6]:
                                    individual = individual + '\t1'
                                else:
                                    individual = individual + '\t0'
                            else:
                                if len(list(df_time_ind_index['allele_snp_id'].unique())) < 2:
                                    if maf_dict[index][1] == df_time_ind_index.iloc[0, 6]:
                                        individual = individual + '\t1'
                                    else:
                                        individual = individual + '\t0'
                                else:
                                    individual = individual + '\t0.5'
                    individual = individual + '\n'
                    g.write(individual)            
                g.close()

            else:
                 sys.exit('Error: provide valid time(s)')
            f.close()
        else:
             sys.exit('Error: provide a minfreq between 0 and 0.5')  
    else:
       sys.exit('Error: provide a valid path to the input file')

if __name__ == '__main__':
     TheRealMcCoilFormat(args.inputfile, args.time, args.minfreq)

