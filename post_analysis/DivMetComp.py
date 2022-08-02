#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:09:26 2022
@author: Frederic Labbe
This script compares the varmodel3 output databases (in SQLite3 format) located in two distinct directories.
It could be used to evaluate new versions of the model, e.g., outputs before vs. after a new implementation.
Four major diversity metrics are compared: the average prevalence, pairwise type sharing (PTS), and numner of strains and genes per replicate.
For each diversity metric, the script plots its distributions and perform a two-sample Kolmogorov-Smirnov test.
It uses the "sampled_hosts", "gene_strain_counts", "sampled_infections", and "sampled_infection_genes" tables from the output database.
Note: A minimum of 2 replicates is required, but we recommend using at least 10 replicates.
usage: python DivMetComp.py --directory1 '/path/to/directory1' --directory2 '/path/to/directory2' --replicates 10
"""

import os.path
import sqlite3
import pandas as pd
import argparse
from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
pd.options.mode.chained_assignment = None
parser = argparse.ArgumentParser()
parser.add_argument('-d1', "--directory1", required = True, help = 'Path to the directory containing the 1st set of varmodel3 output files')
parser.add_argument('-d2', "--directory2", required = True, help = 'Path to the directory containing the 2nd set of varmodel3 output files')
parser.add_argument('-r', "--replicates", type = int, required = True, help = 'Number of replicates')
args = parser.parse_args()

def ModelComparison(directory1, directory2, replicates):
    if replicates > 1:
        ls1 = []; ls2 = []; ls3 = []; ls4 = []
        ls5 = []; ls6 = []; ls7 = []; ls8 = [];
        for replicate in range(1, replicates + 1):
            
            # Extract/calcul average prevalence per replicate
            d1 = []; d2 = []
            df1 = ExtractPrevalence(directory1, replicate)
            df2 = ExtractPrevalence(directory2, replicate)
            if Times(df1, df2):
                for t in df1.time.unique():
                    preval1 = Prevalence(df1, t)
                    preval2 = Prevalence(df2, t)
                    d1.append(preval1)
                    d2.append(preval2)
                ls1.append(np.array(d1).mean())
                ls2.append(np.array(d2).mean())
            else:
                print('Error: the prevalences are not comparable')
                
            # Extract/calcul average gene and strain counts per replicate
            gs1 = ExtractGeneStrains(directory1, replicate)
            gs2 = ExtractGeneStrains(directory2, replicate)
            if Times(gs1, gs2):
                ls3.append(gs1['n_circulating_genes_blood'].mean())
                ls4.append(gs2['n_circulating_genes_blood'].mean())
                ls5.append(gs1['n_circulating_strains_blood'].mean())
                ls6.append(gs2['n_circulating_strains_blood'].mean())
            else:
                print('Error: the prevalences are not comparable')
                
            # Extract/calcul average PTS per replicate
            d3 = []; d4 = []
            inf1 = ExtractPTS(directory1, replicate)
            inf2 = ExtractPTS(directory2, replicate)            
            if Times(inf1, inf2):
                for t in inf1.time.unique():
                    pts1 = PTS(inf1, t)
                    pts2 = PTS(inf2, t)
                    d3.extend(pts1)
                    d4.extend(pts2)
                ls7.append(np.array(d3).mean())
                ls8.append(np.array(d4).mean())
            else:
                print('Error: the prevalences are not comparable')
                
        # Plot histograms and perform two-sample Kolmogorov-Smirnov tests
        f = open('Kolmogorov_Smirnov_test.txt', 'w')
        f.write("measure\tstatistic\tpvalue\n")
        prevalence = [ls1, ls2]; circulating_genes = [ls3, ls4]
        circulating_strains = [ls5, ls6]; pts = [ls7, ls8]
        for measure in ['prevalence', 'circulating_genes', 'circulating_strains', 'pts']:
            plt.hist(eval(measure)[0], alpha = 0.5, color = 'r', label = 'Model 1')
            plt.hist(eval(measure)[1], alpha = 0.5, color = 'g', label = 'Model 2')
            plt.xlabel(measure); plt.ylabel('Count'); plt.legend()
            plt.savefig("{}_hist.pdf".format(measure))
            plt.close()
            kstest = stats.kstest(eval(measure)[0], eval(measure)[1])
            f.write("{}\t{}\t{}\n".format(measure, kstest.statistic, kstest.pvalue))
            if kstest.pvalue < 0.05:
                print('KS test: we reject the null hypothesis for {}'.format(measure))
            else:
                print('KS test: we cannot reject the null hypothesis for {}'.format(measure))
        f.close()
    else:
        print('Error: at least two replicates are required')

def ExtractPrevalence(directory, replicate):
    if os.path.exists("{}/r{}/output.sqlite".format(directory, replicate)):
        con = sqlite3.connect("{}/r{}/output.sqlite".format(directory, replicate))
        df = pd.read_sql_query('SELECT time, id, n_infections_active FROM sampled_hosts', con)
        con.close()
        return df
    else:
       print('Error: provide a valid path to the output files')

def ExtractGeneStrains(directory, replicate):
    if os.path.exists("{}/r{}/output.sqlite".format(directory, replicate)):
        con = sqlite3.connect("{}/r{}/output.sqlite".format(directory, replicate))
        df = pd.read_sql_query('SELECT time, n_circulating_genes_blood, n_circulating_strains_blood FROM gene_strain_counts', con)
        con.close()
        return df
    else:
       print('Error: provide a valid path to the output files')

def ExtractPTS(directory, replicate):
    if os.path.exists("{}/r{}/output.sqlite".format(directory, replicate)):
        con = sqlite3.connect("{}/r{}/output.sqlite".format(directory, replicate))
        df1 = pd.read_sql_query('SELECT time, host_id, infection_id, strain_id, expression_index FROM sampled_infections', con)
        df1.columns = ['time', 'host_id', 'infection_id', 'strain_id', 'expression_index_infection']
        df1 = df1.dropna()
        df2 = pd.read_sql_query('SELECT infection_id, expression_index, allele_id_1, allele_id_2 FROM sampled_infection_genes', con)
        df2.columns = ['infection_id', 'expression_index_gene', 'allele_id_1', 'allele_id_2']
        df = df1.merge(df2, left_on = 'infection_id', right_on = 'infection_id')
        con.close()
        return df
    else:
       print('Error: provide a valid path to the output files')

def Times(df1, df2):
    times1 = df1.time.unique(); times1.sort()
    times2 = df2.time.unique(); times2.sort()
    if np.array_equal(times1, times2):
        return True
    else:
        return False

def Prevalence(df, time):
    part_df = df[df['time'] == time]
    nbpart = len(set(part_df['id']))
    host_df = part_df[part_df['n_infections_active'] > 0]
    nbhost = len(set(host_df['id']))
    preval = nbhost / nbpart
    return preval

def PTS(df, time):    
    # Merge and filter:
    if df['time'].isin([time]).any():
        df_time = df[df['time'] == time]
        df_time["gene_id"] = df_time["allele_id_1"].astype(str) + '_' + df_time["allele_id_2"].astype(str)
        # Convert data between wide and long forms (matrix output; i.e. 0 or 1 for absent or present gene in that strain).
        g = df_time.groupby('infection_id')['gene_id'].apply(list).reset_index()
        genemat_time = g.join(pd.get_dummies(g['gene_id'].apply(pd.Series).stack()).sum(level = 0)).drop('gene_id', 1)
        genemat_time = genemat_time.iloc[: , 1:]
        genemat_time = genemat_time.to_numpy()
        # Cross-product of the transpose of the matrix, divided by the number of genes per strain:
        m0 = (genemat_time > 0).astype(int)
        newmat = m0 @ (m0.T)
        networkSim_time = newmat / (genemat_time > 0).sum(axis = 1)
        networkSim_time_nodiag = networkSim_time[~np.eye(networkSim_time.shape[0], dtype = bool)].reshape(networkSim_time.shape[0], -1)
        networkSim_time_nodiag = networkSim_time_nodiag.flatten()
        return networkSim_time_nodiag

if __name__ == '__main__':
     ModelComparison(args.directory1, args.directory2, args.replicates)
