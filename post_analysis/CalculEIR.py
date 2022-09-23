#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:02:17 2022
@author: flabbe
This script calculates the Entomological Inoculation Rate (EIR).
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
usage: python CalculEIR.py --inputfile '/path/to/file.txt' --time 300 --n_hosts 10000 --summary_period 30 t_year 360
"""

import os.path
import sqlite3
import pandas as pd
import argparse
import sys
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file')
parser.add_argument('-t', "--time", type = int, required = True, help = 'Time to make the calculations')
parser.add_argument('-n', "--n_hosts", type = int, required = True, help = 'Number of hosts used in varmodel3')
parser.add_argument('-s', "--summary_period", type = int, required = True, help = 'How often the summary output was written in varmodel3')
parser.add_argument('-y', "--t_year", type = int, required = True, help = 'Number of time units in a year used in varmodel3')
args = parser.parse_args()

def CalculEIR(inputfile, time, n_hosts, summary_period, t_year):
    if os.path.exists(inputfile):
        os.chdir(os.path.dirname(os.path.abspath(inputfile)))

        # Extract and reformat information from the output database
        con = sqlite3.connect(inputfile)
        df = pd.read_sql_query('SELECT time, n_bites, n_infected_bites FROM summary', con)
        con.close()
        
        # Calcul EIR
        if df['time'].isin([time]).any():
            df_time = df[df['time'] == time]
            eir = df_time.iloc[0]['n_infected_bites'] / n_hosts / summary_period * t_year
            n_bites_host_year = df_time.iloc[0]['n_bites'] / n_hosts / summary_period * t_year
            
            # Export results
            outputfile = inputfile.split("/")[-1].split(".")[0] + "_EIR_" + str(time) + "days.txt"
            f = open(outputfile, 'w')
            f.write("time\teir\tn_bites_host_year\n")
            f.write("{}\t{}\t{}\n".format(time, eir, n_bites_host_year))
            f.close()            
        else:
             sys.exit('Error: provide valid time')        
    else:
       sys.exit('Error: provide a valid path to the input file')

if __name__ == '__main__':
     CalculEIR(args.inputfile, args.time, args.n_hosts, args.summary_period, args.t_year)
     
