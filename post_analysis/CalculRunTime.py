#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 16:51:44 2022
@author: Frederic Labbe
This script calculates and summarize the running times per runs and per replicates.
It uses the "summary" and "meta" tables from the varmodel3 output database(s) (in SQLite3 format).
The results will be stored into two output files: one for the running time per run and another one per replicate.
Note: the two output files will have similar content if there is only one replicate per run.
usage: python CalculRunTime.py --directory '/path/to/directory' --runs 20 --replicates 10
"""

import os.path
import sqlite3
import pandas as pd
import statistics
import argparse
import sys
parser = argparse.ArgumentParser()
parser.add_argument('-d', "--directory", required = True, help = 'Path to the directory containing the input file(s)')
parser.add_argument('-n', "--runs", type = int, required = True, help = 'Number of runs')
parser.add_argument('-r', "--replicates", type = int, required = True, help = 'Number of replicates')
args = parser.parse_args()

def CalculRunTime(directory, runs, replicates):
    if os.path.exists(directory):
        os.chdir(os.path.dirname(os.path.abspath(directory)))
        if runs > 0 and replicates > 0:
        
            # Build dictionaries:
            runs = runs + 1
            replicates = replicates + 1
            meta_rep_dict = {}
            meta_run_dict = {}
            summ_rep_dict = {}
            summ_run_dict = {}
            for run in range(1, runs):
                meta_run_dict[run] = [], []
                summ_run_dict[run] = []
                for replicate in range(1, replicates):
                    if os.path.exists("{}/c{}/r{}/output.sqlite".format(directory, run, replicate)):
                        con = sqlite3.connect("{}/c{}/r{}/output.sqlite".format(directory, run, replicate))
                        meta = pd.read_sql_query('SELECT * FROM meta', con)
                        summary = pd.read_sql_query('SELECT time, exec_time FROM summary', con)
                        con.close()
                        elapsed_time = meta.value[1]
                        summary = summary.iloc[1: , :]
                        summ_rep_dict[run, replicate] = summary
                        time_max = max(summary.time)
                        meta_rep_dict[run, replicate] = elapsed_time, time_max
                        meta_run_dict[run][0].append(elapsed_time)
                        meta_run_dict[run][1].append(time_max)
                        summ_run_dict[run].extend(summary.exec_time.values)
                    else:
                        sys.exit('Error: provide a valid number of runs and/or replicates')
            
            # Running time per run:
            f1 = open("output_runningtime_run.txt", 'w')
            f1.write("run\telapsed_time_avg\telapsed_time_med\telapsed_time_sd\ttime_max_avg\ttime_max_med\ttime_max_sd\texec_time_avg\texec_time_med\texec_time_sd\n")
            for run in range(1, runs):
                elapsed_time_avg = statistics.mean(meta_run_dict[run][0])
                elapsed_time_med = statistics.median(meta_run_dict[run][0])
                time_max_avg = statistics.mean(meta_run_dict[run][1])
                time_max_med = statistics.median(meta_run_dict[run][1])
                exec_time_avg = statistics.mean(summ_run_dict[run])
                exec_time_med = statistics.median(summ_run_dict[run])
                exec_time_sd = statistics.stdev(summ_run_dict[run])
                if replicates > 2:
                    elapsed_time_sd = statistics.stdev(meta_run_dict[run][0])
                    time_max_sd = statistics.stdev(meta_run_dict[run][1])                    
                else:
                    elapsed_time_sd = 0
                    time_max_sd = 0
                f1.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(run, elapsed_time_avg, elapsed_time_med, elapsed_time_sd, time_max_avg, time_max_med, time_max_sd, exec_time_avg, exec_time_med, exec_time_sd))
            f1.close()
            
            # Running time per replicate:
            f2 = open("output_runningtime_rep.txt", 'w')
            f2.write("run\treplicate\telapsed_time\ttime_max\texec_time_avg\texec_time_med\texec_time_sd\n")
            for run in range(1, runs):
                for replicate in range(1, replicates):
                    elapsed_time, time_max = meta_rep_dict[run, replicate]
                    exec_time_avg = statistics.mean(summ_rep_dict[run, replicate].exec_time)
                    exec_time_med = statistics.median(summ_rep_dict[run, replicate].exec_time)
                    exec_time_sd = statistics.stdev(summ_rep_dict[run, replicate].exec_time)
                    f2.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(run, replicate, elapsed_time, time_max, exec_time_avg, exec_time_med, exec_time_sd))
            f2.close()
        else:
            sys.exit('Error: provide a valid number of runs and/or replicates')
    else:
        sys.exit('Error: provide a valid path to the directory containing the input file(s)')
        
if __name__ == '__main__':
     CalculRunTime(args.directory, args.runs, args.replicates)
