#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 16:20:01 2022
@author: Frederic Labbe
This script calculates the targets for the model fitting.
These targets correspond to five diversity metrics: the prevalence, MOI, PTS, and number of strains and genes.
It uses the "sampled_host", "sampled_infections", and "sampled_infection_genes" tables from the varmodel3 output database.
To account for measurement error, calculations can be done with a measurement model using the "--measurement" parameter.
If so, it is required to provide a file describing the distribution of the number of gene per monoclonal infection.
Moreover, it is also required to provide a proportion of host to keep for the calculations (e.g. proportion detected by microscopy).
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
Note: the calculations only take into account the active infections.
usage: python Targets.py --inputfile '/path/to/file.txt' --time 300 --measurement --distribution '/path/to/distribution.txt' --prop 0.57 
"""

import os.path
import sqlite3
import pandas as pd
import argparse
import random
import numpy as np
import sys
import warnings
import statistics
import pickle
from scipy.stats import nbinom
from multiprocessing import Manager
warnings.simplefilter(action='ignore', category=FutureWarning)
pd.options.mode.chained_assignment = None


def create_argparser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--inputfile", required = True, help = 'Path to input file')
    parser.add_argument('-t', "--time", type = int, required = True, help = 'Time at which to calculate the targets')
    parser.add_argument('-m', "--measurement", required = False, action = "store_true", help = 'Wether calculations take into account the measurement error accounting for sub-sampling of var genes')
    parser.add_argument('-p', "--prop", type = float, required = False, help = 'Proportion of positively infected hosts to keep in the calculations; this is to account for the lower detection power of microscopy')
    parser.add_argument('-s', "--supplementaryFileForMOIEst", type = str, required = True, help = 'Path to the file which stores information for MOI estimation')
    parser.add_argument('-a', "--aggregate", type = str, required = True, help = 'How to obtain the MOI distribution at the population level from individual MOI estimates; either pooling the maximum a posteriori MOI estimate of each individual host, or using the technique called mixture distribution, i.e., a weighted summation of all individual MOI distributions; pool or mixtureDist')
    return parser


def CalculTargetsMeasFunc(inputfile, time, prop, supplementaryFileForMOIEst, aggregate, lock):
    if os.path.exists(inputfile):
        with lock:
            con = sqlite3.connect(inputfile)
            df1 = pd.read_sql_query('SELECT time, id, n_infections_active, birth_time FROM sampled_hosts', con)
            df2 = pd.read_sql_query('SELECT time, host_id, infection_id, strain_id, expression_index FROM sampled_infections', con)
            df3 = pd.read_sql_query('SELECT infection_id, expression_index, group_id, allele_id_1, allele_id_2 FROM sampled_infection_genes', con)
            con.close()

        df1.columns = ['time', 'host_id', 'n_infections_active', 'birth_time']
        df2.columns = ['time', 'host_id', 'infection_id', 'strain_id', 'expression_index_infection']
        df2 = df2.dropna()
        df3.columns = ['infection_id', 'expression_index_gene', 'group_id', 'allele_id_1', 'allele_id_2']

        if len(df2) > 0:
            # Filter/reformat data
            df1 = Times(df1, time)
            df2 = Times(df2, time)
            df = df2.merge(df3, left_on = 'infection_id', right_on = 'infection_id')
            df["gene_id"] = df["allele_id_1"].astype(str) + '_' + df["allele_id_2"].astype(str)

            # Host subsampling
            sub = round(prop * len(df['host_id'].unique()))
            meas = random.sample(list(df['host_id'].unique()), sub)
            rslt_df = df[df['host_id'].isin(meas)]
            
            # load measurement error
            with open(supplementaryFileForMOIEst,'rb') as handle:
                MOIEstObjs = pickle.load(handle)
            errA = MOIEstObjs[0]
            errBC = MOIEstObjs[1]
            priors = MOIEstObjs[2]
            err_distribution = MOIEstObjs[3]
            probSizeGivenMOIDict = MOIEstObjs[4]
            maxMOI = MOIEstObjs[5]
            repSizeLow = errBC[0].min()
            repSizeHigh = errBC[0].max()
            
            # Calculate the targets
            preval = Prevalence(rslt_df, df1)
            nbstrain = len(set(rslt_df['strain_id']))
            MOIs = pd.DataFrame(columns = ['MOI', 'Prob', 'host_id'])
            genesA = []
            genesBC = []
            genes = []
            subsetsA = pd.DataFrame(columns = ['gene_id', 'strain_id', 'host_id'])
            subsetsBC = pd.DataFrame(columns = ['gene_id', 'strain_id', 'host_id'])
            subsets = pd.DataFrame(columns = ['gene_id', 'strain_id', 'host_id'])

            for host in rslt_df['host_id'].unique():
                var_host = rslt_df[rslt_df['host_id'] == host]
                size = max(var_host['expression_index_gene'])
                subsampA = []
                subsampBC = []
                for strain in var_host['strain_id'].unique():
                    var_strain = var_host[var_host['strain_id'] == strain]
                    
                    varA_strain = var_strain[var_strain['group_id'] == 1]
                    varBC_strain = var_strain[var_strain['group_id'] == 2]
                    
                    nb_varA_strain = len(varA_strain['gene_id'].unique())
                    nb_varBC_strain = len(varBC_strain['gene_id'].unique())
                    
                    sampA = int(random.choices(population = np.array(errA.iloc[:, 0]), weights = np.array(errA.iloc[:, 1]), k = 1)[0])
                    sampBC = int(random.choices(population = np.array(errBC.iloc[:, 0]), weights = np.array(errBC.iloc[:, 1]), k = 1)[0])
                    while (sampA > nb_varA_strain):
                        sampA = int(random.choices(population = np.array(errA.iloc[:, 0]), weights = np.array(errA.iloc[:, 1]), k = 1)[0])
                    while (sampBC > nb_varBC_strain):
                        sampBC = int(random.choices(population = np.array(errBC.iloc[:, 0]), weights = np.array(errBC.iloc[:, 1]), k = 1)[0])
                    varA_samp = random.sample(list(varA_strain.gene_id), sampA)
                    varBC_samp = random.sample(list(varBC_strain.gene_id), sampBC)
                    
                    subsampA.extend(varA_samp)
                    subsampBC.extend(varBC_samp)
                    genes.extend(varA_samp)
                    genes.extend(varBC_samp)
                    genesA.extend(varA_samp)
                    genesBC.extend(varBC_samp)
                    subsetA = pd.DataFrame(varA_samp, columns = ['gene_id'])
                    subsetA['strain_id'] = np.repeat(varA_strain.strain_id.unique(), len(subsetA))
                    subsetA['host_id'] = np.repeat(host, len(subsetA))
                    subsetBC = pd.DataFrame(varBC_samp, columns = ['gene_id'])
                    subsetBC['strain_id'] = np.repeat(varBC_strain.strain_id.unique(), len(subsetBC))
                    subsetBC['host_id'] = np.repeat(host, len(subsetBC))
                    subsets = subsets.append(subsetA)
                    subsets = subsets.append(subsetBC)
                    subsetsA = subsetsA.append(subsetA)
                    subsetsBC = subsetsBC.append(subsetBC)    
                
                nb_varBC_err = len(set(subsampBC))
                MOIInd = MOIest(nb_var = nb_varBC_err, priors = priors, probSizeGivenMOI = probSizeGivenMOIDict, maxMOI = maxMOI, repSizeLow = repSizeLow, repSizeHigh = repSizeHigh)
                MOIInd['host_id'] = np.repeat(host, len(MOIInd))
                MOIs = MOIs.append(MOIInd)
            if aggregate == "pool":
                MOIsPop = MOIs.sort_values(['Prob'],ascending=False).groupby('host_id').head(1)
                MOIvar = np.asarray(MOIsPop['MOI']).mean()
            elif aggregate == "mixtureDist":
                MOIsPop = MOIs.groupby('MOI').agg({'Prob': 'sum'})
                # print(MOIsPop)
                # print(np.dot(np.asarray(MOIsPop.index), np.asarray(MOIsPop['Prob'])))
                # print(sum(np.asarray(MOIsPop['Prob'])))
                MOIvar = np.dot(np.asarray(MOIsPop.index), np.asarray(MOIsPop['Prob']))/sum(np.asarray(MOIsPop['Prob']))

            pts = PTS(subsets)
            ptsA = PTS(subsetsA)
            ptsBC = PTS(subsetsBC)
            nbgene = len(np.unique(genes))
            nbgeneA = len(np.unique(genesA))
            nbgeneBC = len(np.unique(genesBC))
            return preval, MOIvar, pts, ptsA, ptsBC, nbstrain, nbgene, nbgeneA, nbgeneBC
        else:
            return 0, 0, 0, 0, 0, 0, 0, 0, 0
    else:
        sys.exit('Error: provide a valid path to the input file') 
   
def CalculTargetsMeas(inputfile, time, prop, supplementaryFileForMOIEst, aggregate):
    outputfile = inputfile.split("/")[-1].split(".")[0] + "_targets_" + str(time) + "days_meas.txt"
    f = open(outputfile, 'w')
    if prop is not None and (0 < prop <= 1):
        if os.path.exists(supplementaryFileForMOIEst):
            f.write("prevalence_err\tMOI_err\tPTS_err\tPTS_err_A\tPTS_err_BC\tn_strains_err\tn_genes_err\tn_genes_err_A\tn_genes_err_BC\n")
        else:
            sys.exit('Error: provide a valid path to the measurement model')
        with Manager() as manager:
            lock = manager.Lock()
            preval, MOIvar, pts, ptsA, ptsBC, nbstrain, nbgene, nbgeneA, nbgeneBC = CalculTargetsMeasFunc(inputfile, time, prop, supplementaryFileForMOIEst, aggregate, lock)
            f.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(preval, MOIvar, pts, ptsA, ptsBC, nbstrain, nbgene, nbgeneA, nbgeneBC))
            f.close()
    else:
        sys.exit('Error: provide a valid proportion of host to keep in the calculations')


def CalculTargets(inputfile, time):
    outputfile = inputfile.split("/")[-1].split(".")[0] + "_targets_" + str(time) + "days.txt"
    f = open(outputfile, 'w')
    f.write("prevalence\tMOI\tPTS\tn_strains\tn_genes\n")
    if os.path.exists(inputfile):

        # Extract data
        con = sqlite3.connect(inputfile)
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
        
        # Export targets        
        f.write("{}\t{}\t{}\t{}\t{}\n".format(preval, MOIvar, pts, nbstrain, nbgene))
        f.close()
    else:
       sys.exit('Error: provide a valid path to the input file')
 
def Times(df, time):
    if df['time'].isin([time]).any():
        df = df[df['time'] == time]
        return df
    else:
        sys.exit('Error: provide a valid time')

def Prevalence(df1, df2):
    nbpart = len(set(df2['host_id']))
    nbhost = len(set(df1['host_id']))
    preval = nbhost / nbpart
    return preval
   
def PTS(df):
    # if len(df['host_id'].unique()) > 1000:
    #     select = np.random.choice(df['host_id'].unique(), size = 1000)
    #     df = df[df['host_id'].isin(select)] 
    g = df.groupby('host_id')['gene_id'].apply(list).reset_index()
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
        
def MOIest(nb_var, priors, probSizeGivenMOI, maxMOI = 20, repSizeLow = 10, repSizeHigh = 45):
    if nb_var < repSizeLow:
        warnings.warn("The number of non-upsA DBLa type is fewer than the lower number presented in your repertoire size distribution. Automatically assign 1 to be the MOI value.")
        data_temp = {"MOI": list(range(1, maxMOI + 1, 1)), "Prob": [1.0] + [0.0]*(maxMOI - 1)}
        probMOIGivenIsoSize = pd.DataFrame(data_temp)
    elif nb_var > repSizeHigh*maxMOI:
        warnings.warn("The number of non-upsA DBLa type is greater than the maximum possible one, i.e., the maximum repertoire size * maxMOI. Automatically assign maxMOI to be the MOI value.")
        data_temp = {"MOI": list(range(1, maxMOI + 1, 1)), "Prob": [0.0]*(maxMOI - 1) + [1.0]}
        probMOIGivenIsoSize = pd.DataFrame(data_temp)
    elif nb_var >= repSizeLow and nb_var <= repSizeHigh*maxMOI:
        denominator = 0.0
        numerators = []
        for MOI in range(1, maxMOI+1, 1): 
            prob1_temp = probSizeGivenMOI[MOI]
            if any(prob1_temp["Size"] == nb_var):
                prob1 = prob1_temp["Prob"][prob1_temp["Size"] == nb_var].tolist()[0]
            else:
                prob1 = 0.0
            prob2 = priors["Prob"][priors["MOI"] == MOI].tolist()[0]
            denominator_temp = prob1 * prob2
            numerators.append(denominator_temp)
            denominator += denominator_temp
        if denominator == 0.0:
            sys.exit("Ah! The denominator is zero! No MOI values can return the isolate size observed!")
        probMOIGivenIsoSize = [x/denominator for x in numerators]
    MOI_temp = pd.DataFrame({'MOI': list(range(1, maxMOI + 1, 1)), 'Prob': probMOIGivenIsoSize})
    return MOI_temp

# if __name__ == '__main__':
#     parser = create_argparser()
#     args = parser.parse_args()
#     if args.measurement:
#         CalculTargetsMeas(args.inputfile, args.time, args.prop, args.supplementaryFileForMOIEst, args.aggregate)
#         # CalculTargetsMeas(args.inputfile, args.time, args.distribution, args.prop)
#     else:
#         CalculTargets(args.inputfile, args.time)
