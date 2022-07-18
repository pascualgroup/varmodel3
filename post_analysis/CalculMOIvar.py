#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 15:10:13 2022
@author: Frederic Labbe
This script calculates the multiplicity of infection (MOI) per hosts using the varcoding approach.
For each host, it also exports the number of unique var genes and the true MOI.
It uses the "sampled_infections", "sample_hosts", and "sampled_infection_genes" tables from the varmodel3 output database (in SQLite3 format).
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
Note: the MOI calculations only take into account the active infection(s).
usage: python CalculMOIvar.py --inputfile '/path/to/file.txt' --time 300
"""

import os.path
import sqlite3
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file')
parser.add_argument('-t', "--time", type = int, required = True, help = 'Time to make the calculations')
args = parser.parse_args()

host = 4315
def CalculMOIvar(inputfile, time):
    if os.path.exists(inputfile):
        os.chdir(os.path.dirname(os.path.abspath(inputfile)))
        
        # Extract and reformat information from the output database:
        con = sqlite3.connect('output.sqlite')
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
            
            # Calculate the MOI using the varcoding approach:
            df_time["gene_id"] = df_time["allele_id_1"].astype(str) + '_' + df_time["allele_id_2"].astype(str)
            outputfile = inputfile.split("/")[-1].split(".")[0] + "_MOIvar_" + str(time) + "days.txt"
            f = open(outputfile, 'w')
            f.write("host_id\tnb_var\tMOIvar\tMOI\n")
            for host in df_time['host_id'].unique():
                var = df_time[df_time['host_id'] == host]
                size = max(var['expression_index_gene'])
                nb_var = len(var['gene_id'].unique())
                temp = nb_var
                MOI = int(df3[df3['host_id'] == host].n_infections_active)
                MOIvar = 1
                while (temp > size):
                    temp = temp - size
                    MOIvar += 1
                f.write("{}\t{}\t{}\t{}\n".format(host, nb_var, MOIvar, MOI))        
            f.close()
        else:
             print('Error: provide a valid time')
    else:
       print('Error: provide a valid path to the input file') 

if __name__ == '__main__':
     CalculMOIvar(args.inputfile, args.time)
