#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 14:23:08 2023
@author: Frederic Labbe
This script calculates the gene frequency distribution.
It uses the "sampled_infections"  and "sampled_infection_genes" tables from the varmodel3 output database (in SQLite3 format).
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
Note: the calculations only take into account the active infection(s).
usage: python GeneFrequencyDistribution.py --inputfile '/path/to/file.txt' --time 300
"""

import os.path
import sqlite3
import pandas as pd
import sys
import argparse
pd.options.mode.chained_assignment = None
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file')
parser.add_argument('-t', "--time", type = int, required = True, help = 'Time to make the calculations')
args = parser.parse_args()

def GeneFrequencyDistribution(inputfile, time):
    if os.path.exists(inputfile):
        os.chdir(os.path.dirname(os.path.abspath(inputfile)))

        # Extract and reformat information from the output database
        con = sqlite3.connect(inputfile)
        df1 = pd.read_sql_query('SELECT time, host_id, infection_id, strain_id, expression_index FROM sampled_infections', con)
        df1.columns = ['time', 'host_id', 'infection_id', 'strain_id', 'expression_index_infection']
        df1 = df1.dropna()
        df2 = pd.read_sql_query('SELECT infection_id, expression_index, allele_id_1, allele_id_2 FROM sampled_infection_genes', con)
        df2.columns = ['infection_id', 'expression_index_gene', 'allele_id_1', 'allele_id_2']
        con.close()
        
        # Merge and filter
        df = df1.merge(df2, left_on = 'infection_id', right_on = 'infection_id')
        if df['time'].isin([time]).any():
            df_time = df[df['time'] == time]
            df_time["gene_id"] = df_time["allele_id_1"].astype(str) + '_' + df_time["allele_id_2"].astype(str)
            
            # Calculate the gene frequency distribution
            frequency = df_time["gene_id"] .value_counts().rename_axis('gene_id').reset_index(name = 'frequency')
            outputfile = inputfile.split("/")[-1].split(".")[0] + "_GenfreqDist_" + str(time) + "days.csv"
            frequency.to_csv(outputfile, index = False, header = True)    
        else:
             sys.exit('Error: provide a valid time')
    else:
       sys.exit('Error: provide a valid path to the input file')
    
if __name__ == '__main__':
     GeneFrequencyDistribution(args.inputfile, args.time)
