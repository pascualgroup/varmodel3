#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 15:41:27 2023
@author: Frederic Labbe
This script calculates the Shannon and Simpson diversity indexes of the var genes.
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
usage: python CalculDivIndex.py --inputfile '/path/to/file.txt' --time 300 --prop 0.57
"""

import os.path
import sqlite3
import pandas as pd
import argparse
import random
import sys
from math import log as ln
pd.options.mode.chained_assignment = None
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file')
parser.add_argument('-t', "--time", type = int, required = True, help = 'Time to make the calculations')
parser.add_argument('-p', "--prop", type = float, required = False, help = 'Proportion of host to keep in the calculations')
args = parser.parse_args()

def CalculDivIndex(inputfile, time, prop):
    if os.path.exists(inputfile):
        os.chdir(os.path.dirname(os.path.abspath(inputfile)))
        
        # Extract and reformat information from the output database:
        con = sqlite3.connect(inputfile)
        df1 = pd.read_sql_query('SELECT time, host_id, infection_id, strain_id, expression_index FROM sampled_infections', con)
        df1.columns = ['time', 'host_id', 'infection_id', 'strain_id', 'expression_index_infection']
        df1 = df1.dropna()
        df2 = pd.read_sql_query('SELECT infection_id, expression_index, allele_id_1, allele_id_2 FROM sampled_infection_genes', con)
        df2.columns = ['infection_id', 'expression_index_gene', 'allele_id_1', 'allele_id_2']
        df3 = pd.read_sql_query('SELECT time, id, n_infections_active, birth_time FROM sampled_hosts', con)
        df3.columns = ['time', 'host_id', 'n_infections_active', 'birth_time']
        con.close()

        # Merge and filter:
        df = df1.merge(df2, left_on = 'infection_id', right_on = 'infection_id')
        if df['time'].isin([time]).any():
            df_time = df[df['time'] == time]
            df3 = df3[df3['time'] == time]
            
            # Host subsampling:
            df_time["gene_id"] = df_time["allele_id_1"].astype(str) + '_' + df_time["allele_id_2"].astype(str)
            if prop is not None:
                if (0 < prop <= 1):
                    sub = round(prop * len(df_time['host_id'].unique()))
                    meas = random.sample(list(df_time['host_id'].unique()), sub)
                    rslt_df_time = df_time[df_time['host_id'].isin(meas)]
                else:
                    sys.exit('Error: provide a valid proportion of host to keep in the calculations')
            
            # Calculate the Shannon diversity index:
            def Shannon(n, N):
                """ Relative abundance """
                if n == 0:
                    return 0
                else:
                    return (float(n)/N) * ln(float(n)/N)
            if prop is not None and (0 < prop <= 1):
                counts_err = rslt_df_time['gene_id'].value_counts().tolist()
                N_err = rslt_df_time["gene_id"].value_counts().sum()
                sdi_err = -sum(Shannon(n, N_err) for n in counts_err if n != 0)
            counts = df_time['gene_id'].value_counts().tolist()
            N = df_time["gene_id"].value_counts().sum()            
            sdi = -sum(Shannon(n, N) for n in counts if n != 0)
            
            # Calculate the Simpson diversity and inverse Simpson diversity indexes:
            def Simpson(n, N):
                """ Relative abundance """
                if n == 0:
                    return 0
                else:
                    return float(n)/N
            if prop is not None and (0 < prop <= 1):
                sidi_err = sum(Simpson(n, N_err)**2 for n in counts_err if n != 0)
                sidi_err_inv = float(1)/sidi_err 
            sidi = sum(Simpson(n, N)**2 for n in counts if n != 0)
            sidi_inv = float(1)/sidi
            
            # Export results:
            outputfile = inputfile.split("/")[-1].split(".")[0] + "_SDI_" + str(time) + "days.txt"
            f = open(outputfile, 'w')
            if prop is not None and (0 < prop <= 1):
                f.write("Shannon_Diversity_Index\tShannon_Diversity_Index_Err\tSimpson_Diversity_Index\tSimpson_Diversity_Index_Err\tInverse_Simpson_Diversity_Index\tInverse_Simpson_Diversity_Index_Err\n")
                f.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sdi, sdi_err, sidi, sidi_err, sidi_inv, sidi_err_inv))
            else:
                f.write("Shannon_Diversity_Index\tSimpson_Diversity_Index\tInverse_Simpson_Diversity_Index\n")
                f.write("{}\t{}\t{}\n".format(sdi, sidi, sidi_inv))
            f.close()
        else:
             sys.exit('Error: provide a valid time')
    else:
       sys.exit('Error: provide a valid path to the input file')

if __name__ == '__main__':
     CalculDivIndex(args.inputfile, args.time, args.prop)
