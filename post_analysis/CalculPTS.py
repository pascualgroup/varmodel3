#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:10:13 2022
@author: Frederic Labbe
This script calculates the pairwise type sharing (PTS) between var repertoires.
PTSij = 2nij / (ni + nj), where ni and nj are the number of unique var genes within each repertoire i and j,
and nij is the total number of var genes shared between repertoires i and j.
It uses the "sampled_infections", and "sampled_infection_genes" tables from the varmodel3 output database (in SQLite3 format).
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
Note: the MOI calculations only take into account the active infection(s).
usage: python CalculPTS.py --inputfile '/path/to/file.txt' --time 300
"""

import os.path
import sqlite3
import pandas as pd
import argparse
import sys
pd.options.mode.chained_assignment = None
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file')
parser.add_argument('-t', "--time", type = int, required = True, help = 'Time to make the calculations')
args = parser.parse_args()

def CalculPTS(inputfile, time):
    if os.path.exists(inputfile):
        os.chdir(os.path.dirname(os.path.abspath(inputfile)))
        
        # Extract and reformat information from the output database:
        con = sqlite3.connect(inputfile)
        df1 = pd.read_sql_query('SELECT time, host_id, infection_id, strain_id, expression_index FROM sampled_infections', con)
        df1.columns = ['time', 'host_id', 'infection_id', 'strain_id', 'expression_index_infection']
        df1 = df1.dropna()
        df2 = pd.read_sql_query('SELECT infection_id, expression_index, allele_id_1, allele_id_2 FROM sampled_infection_genes', con)
        df2.columns = ['infection_id', 'expression_index_gene', 'allele_id_1', 'allele_id_2']
        df = df1.merge(df2, left_on = 'infection_id', right_on = 'infection_id')
        con.close()
        
        # Merge and filter:
        if df['time'].isin([time]).any():
            df_time = df[df['time'] == time]
            df_time["gene_id"] = df_time["allele_id_1"].astype(str) + '_' + df_time["allele_id_2"].astype(str)

            # Convert data between wide and long forms (matrix output; i.e. 0 or 1 for absent or present gene in that strain).
            g = df_time.groupby('strain_id')['gene_id'].apply(list).reset_index()
            genemat_time = g.join(pd.get_dummies(g['gene_id'].apply(pd.Series).stack()).sum(level = 0)).drop('gene_id', 1)
            strain_id = genemat_time['strain_id']
            genemat_time = genemat_time.iloc[: , 1:]
            genemat_time = genemat_time.to_numpy()

            # Cross-product of the transpose of the matrix, divided by the number of genes per strain:
            m0 = (genemat_time > 0).astype(int)
            newmat = m0 @ (m0.T)
            networkSim_time = newmat / (genemat_time > 0).sum(axis = 1)

            # Export the results:
            outputfile = inputfile.split("/")[-1].split(".")[0] + "_PTS_" + str(time) + "days.csv"
            out = pd.DataFrame(networkSim_time)
            out = out.rename(columns = strain_id)
            out.index = strain_id
            out.to_csv(outputfile, index = True, header = True)
        
        else:
            sys.exit('Error: provide a valid time')
    else:
       sys.exit('Error: provide a valid path to the input file')

if __name__ == '__main__':
     CalculPTS(args.inputfile, args.time)
