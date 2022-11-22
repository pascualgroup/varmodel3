#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:34:30 2022
@author: Frederic Labbe
This script calculates the targets for the model fitting.
Four major diversity metrics are compared: the prevalence, multiplicity of infection (MOI), pairwise type sharing (PTS), and number of strains and genes.
It uses the "sampled_host", "sampled_infections", and "sampled_infection_genes" tables from the varmodel3 output database (in SQLite3 format).
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
Note: the calculations only take into account the active infections.
usage: python Targets.py --inputfile '/path/to/file.txt' --time 300 --measurement '/path/to/measurement.txt'
"""

import os.path
import sqlite3
import pandas as pd
import argparse
import random
import numpy as np
pd.options.mode.chained_assignment = None


def create_args_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file')
    parser.add_argument('-t', "--time", type = int, required = True, help = 'Time to make the calculate the targets')
    parser.add_argument('-m', "--measurement", required = False, help = 'Path to the measurement model')
    return parser


def calc_targets(input_file, time, measurement):
    # Extract data
    con = sqlite3.connect(input_file)
    df1 = pd.read_sql_query('SELECT time, id, n_infections_active, birth_time FROM sampled_hosts', con)
    df1.columns = ['time', 'host_id', 'n_infections_active', 'birth_time']
    df2 = pd.read_sql_query('SELECT time, host_id, infection_id, strain_id, expression_index FROM sampled_infections', con)
    df2.columns = ['time', 'host_id', 'infection_id', 'strain_id', 'expression_index_infection']
    df2 = df2.dropna()
    df3 = pd.read_sql_query('SELECT infection_id, expression_index, allele_id_1, allele_id_2 FROM sampled_infection_genes', con)
    df3.columns = ['infection_id', 'expression_index_gene', 'allele_id_1', 'allele_id_2']
    con.close()
    
    # Filter/reformat data
    df1 = Times(df1, time)
    df2 = Times(df2, time)
    df = df2.merge(df3, left_on = 'infection_id', right_on = 'infection_id')
    df["gene_id"] = df["allele_id_1"].astype(str) + '_' + df["allele_id_2"].astype(str)

    if measurement is not None:
        # Host subsampling
        sub = round(0.57 * len(df['host_id'].unique()))
        meas = random.sample(list(df['host_id'].unique()), sub)
        rslt_df = df[df['host_id'].isin(meas)]

        # Calculate the targets
        preval = Prevalence(rslt_df, df1)
        nbstrain = len(set(rslt_df['strain_id']))
        MOIs = []
        genes = []
        dfs = [] # subsets = pd.DataFrame(columns = ['strain_id', 'gene_id'])
        for host in rslt_df['host_id'].unique():
            var = rslt_df[rslt_df['host_id'] == host]
            size = max(var['expression_index_gene'])
            err = pd.read_csv(measurement, sep = '\t').values.tolist()
            err = pd.DataFrame(err)
            subsamp = []
            for strain in var['strain_id'].unique():
                var_strain = var[var['strain_id'] == strain]
                nb_var_strain = len(var_strain['gene_id'].unique())
                samp = int(random.choices(population = np.array(err.iloc[:, 0]), weights = np.array(err.iloc[:, 1]), k = 1)[0])
                while (samp > nb_var_strain):
                        samp = int(random.choices(population = np.array(err.iloc[:, 0]), weights = np.array(err.iloc[:, 1]), k = 1)[0])
                var_samp = random.sample(list(var_strain.gene_id), samp)
                subsamp.extend(var_samp)
                genes.extend(var_samp)
                subset = pd.DataFrame(var_samp, columns = ['gene_id'])
                subset['strain_id'] = np.repeat(var_strain.strain_id.unique(), len(subset))            
                dfs.append(subset)
            nb_var_err = len(set(subsamp))
            MOI = 1
            while (nb_var_err > size):
                nb_var_err = nb_var_err - size
                MOI += 1
            MOIs.append(MOI)
        MOIvar = np.asarray(MOIs).mean()
        subsets = pd.concat(dfs)
        pts = PTS(subsets)
        nbgene = len(np.unique(genes))
        
    # Otherwise        
    else:
        # Calculate the targets
        preval = Prevalence(df, df1)
        nbstrain = len(set(df['strain_id']))
        MOIs = []
        genes = []
        for host in df['host_id'].unique():
            var = df[df['host_id'] == host]
            size = max(var['expression_index_gene'])
            nb_var = len(var['gene_id'].unique())
            var_list = var['gene_id'].unique()
            genes.extend(var_list)
            MOI = 1
            while (nb_var > size):
                nb_var = nb_var - size
                MOI += 1
            MOIs.append(MOI)
        MOIvar = np.asarray(MOIs).mean()
        pts = PTS(df)
        nbgene = len(np.unique(genes))
    

    return (preval, MOIvar, pts, nbstrain, nbgene)

def CalculTargets(inputfile, time, measurement):
    outputfile = inputfile.split("/")[-1].split(".")[0] + "_targets_" + str(time) + "days.txt"
    f = open(outputfile, 'w')
    if args.measurement is not None:
        if os.path.exists(measurement):
            f.write("prevalence_err\tMOI_err\tPTS_err\tn_strains_err\tn_genes_err\n")
        else:
            print('Error: provide a valid path to the measurement model') 
    else:
        f.write("prevalence\tMOI\tPTS\tn_strains\tn_genes\n")    

    if os.path.exists(inputfile):
        preval, MOIvar, pts, nbstrain, nbgene = calc_targets(inputfile, time, measurement)
        # Export targets        
        f.write("{}\t{}\t{}\t{}\t{}\n".format(preval, MOIvar, pts, nbstrain, nbgene))
        f.close()        
    else:
       print('Error: provide a valid path to the input file')

def Times(df, time):
    if df['time'].isin([time]).any():
        df = df[df['time'] == time]
        return df
    else:
        print('Error: provide a valid time')

def Prevalence(df1, df2):
    nbpart = len(set(df2['host_id']))
    nbhost = len(set(df1['host_id']))
    preval = nbhost / nbpart
    return preval
   
def PTS(df):
    if len(df['strain_id'].unique()) > 1000:
        df = df.sample(n = 1000)
    g = df.groupby('strain_id')['gene_id'].apply(list).reset_index()
    genemat = g.join(pd.get_dummies(g['gene_id'].apply(pd.Series).stack()).sum(level = 0)).drop('gene_id', 1)
    genemat = genemat.iloc[: , 1:]
    genemat = genemat.to_numpy()
    m0 = (genemat > 0).astype(int)
    newmat = m0 @ (m0.T)
    networkSim = newmat / (genemat > 0).sum(axis = 1)
    networkSim_nodiag = networkSim[~np.eye(networkSim.shape[0], dtype = bool)].reshape(networkSim.shape[0], -1)
    networkSim_nodiag = networkSim_nodiag.flatten()
    pts = networkSim_nodiag.mean()
    return pts

if __name__ == '__main__':
    parser = create_args_parser()
    args = parser.parse_args()
    CalculTargets(args.inputfile, args.time, args.measurement)
