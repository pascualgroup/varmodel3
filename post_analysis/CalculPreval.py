#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 11 16:34:30 2022
@author: Frederic Labbe
This script calculates the prevalence, number of participants, and number of hosts.
It uses the "sampled_host" table from the varmodel3 output database (in SQLite3 format).
Note: the input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
Note: the number of hosts and prevalence only take into account the hosts with at least one active infection.
usage: python CalculPreval.py --inputfile '/path/to/file.txt' --time 30
"""

import os.path
import sqlite3
import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-i', "--inputfile", required = True, help = 'Path to the input file')
parser.add_argument('-t', "--time", type = int, required = True, help = 'Time to make the calculations')
args = parser.parse_args()

def CalculPreval(inputfile, time):
    if os.path.exists(inputfile):
        os.chdir(os.path.dirname(os.path.abspath(inputfile)))
        con = sqlite3.connect(inputfile)
        df = pd.read_sql_query('SELECT * FROM sampled_hosts', con)
        con.close()
        if df['time'].isin([time]).any():
            part_df = df[df['time'] == time]
            host_df = part_df[part_df['n_infections_active'] > 0]
            nbpart = len(set(part_df['id']))
            nbhost = len(set(host_df['id']))
            preval = nbhost / nbpart
            outputfile = inputfile.split("/")[-1].split(".")[0] + "_prevalence_" + str(time) + "days.txt"
            f = open(outputfile, 'w')
            f.write("time\tnb_participants\tnb_hosts\tprevalence\n")
            f.write("{}\t{}\t{}\t{}\n".format(time, nbpart, nbhost, preval))        
            f.close()
        else:
            print('Error: provide a valid time')
    else:
       print('Error: provide a valid path to the input file')     

if __name__ == '__main__':
     CalculPreval(args.inputfile, args.time)
