import pandas as pd
import sqlite3
import os


def run(input_file):
    con = sqlite3.connect(input_file)
    df1 = pd.read_sql_query('SELECT time, id, n_infections_active, birth_time FROM sampled_hosts', con)
    df2 = pd.read_sql_query('SELECT time, host_id, infection_id, strain_id, expression_index FROM sampled_infections', con)
    df3 = pd.read_sql_query('SELECT infection_id, expression_index, group_id, allele_id_1, allele_id_2 FROM sampled_infection_genes', con)
    con.close()

    df1.columns = ['time', 'host_id', 'n_infections_active', 'birth_time']
    df2.columns = ['time', 'host_id', 'infection_id', 'strain_id', 'expression_index_infection']
    df2 = df2.dropna()
    df3.columns = ['infection_id', 'expression_index_gene', 'group_id', 'allele_id_1', 'allele_id_2']

    parent = os.path.dirname(input_file)
    df1_f = os.path.join(parent, 'df1.pq')
    df2_f = os.path.join(parent, 'df2.pq')
    df3_f = os.path.join(parent, 'df3.pq')

    df1.to_parquet(df1_f)
    df2.to_parquet(df2_f)
    df3.to_parquet(df3_f)

    return [df1_f, df2_f, df3_f]
