# Post analysis

This is a collection of scripts to analyze the sqlite output. Below are notes about some useful tools. Not everything is documented yet, but most scripts have some helpful information if you type `python script.py -h`.

## Contents

* [Calculate MOIvar](#Calculate-MOIvar)
* [Calculate PTS](#Calculate-PTS)
* [Calculate running times](#Calculate-running-times)
* [Calculate prevalence](#Calculate-prevalence)
* [Calculate SNP call proportions](#Calculate-SNP-call-proportions)
* [Calculate targets](#Calculate-targets)
* [Compare diversity metrics](#Compare-diversity-metrics)
* [Convert SNPs in THE REAL McCOIL format](#Convert-SNPs-in-THE-REAL-McCOIL-format)
* [SNP measurement model](#SNP-measurement-model)

## Calculate MOIvar
This script calculates the multiplicity of infection (MOI) per hosts using the *var*coding approach. The *var*coding approach (also termed *var* genotyping or *var* fingerprinting), employs the highly polymorphic sequences encoding the immunogenic DBLα domain of PfEMP1 (*Plasmodium falciparum* erythrocyte membrane protein 1), the major surface antigen of the blood stage of infection ([Rask *et al.* 2010](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000933)). The multigene family known as *var* encodes variants of this surface antigen which can reach tens of thousands of variants in endemic populations ([Tonkin-Hill *et al.* 2021](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009269)). The extensive diversity of the *var* gene family together with the very low percentage of *var* genes shared between parasites facilitate measuring MOI by amplifying, pooling, sequencing, and counting the number of DBLα types in a host.

For each host, this script also exports from the simulations the number of unique *var* genes and the true MOI. To account for *var* gene potential sampling errors, an optional measurement model is also available, which sub-sampled the number of *var* genes per strains according to a weight reflecting the *var* gene counts density function. These weights are provided in a file which must contain two columns: 1) the number of *var* genes (e.g. 40, ..., 45), and 2) the weights (e.g. 0.001, ..., 0.075). The script uses the `sampled_infections`, `sample_hosts`, and `sampled_infection_genes` tables from the varmodel3 output database (in [SQLite3 format](https://www.sqlite.org/fileformat.html)).

#### Example command
`python CalculMOIvar.py --inputfile '/path/to/file.txt' --time 300 --measurement '/path/to/measurement.txt'`

`python CalculMOIvar.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `inputfile` | Path to the input file (required) |
| `time`  | Time to make the calculations (required) |
| `measurement`  | Path to the the measurement model (optional) |
#### Notes
The MOI calculations only take into account the active infection(s). The input file name should not contain a `.` except before the extension (*e.g.,* `input_file_name.sqlite`). To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [random](https://docs.python.org/3/library/random.html), [numpy](https://numpy.org/), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

## Calculate PTS
This script calculates the pairwise type sharing (PTS) between var repertoires. *PTSij = 2nij / (ni + nj)*, where *ni* and *nj* are the number of unique *var* genes within each repertoire *i* and *j*, and *nij* is the total number of *var* genes shared between repertoires *i* and *j* [(Barry *et al.* 2007)](https://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.0030034). When the `varmodel3` simulations involved very high diversity (*e.g.,* `n_genes_initial > 20000)`, the PTS is only calculated for a subsample of 1,000 strains to avoid potential memory issues. It uses the `sampled_infections`, and `sampled_infection_genes` tables from the varmodel3 output database (in [SQLite3 format](https://www.sqlite.org/fileformat.html)).
#### Example command
`python CalculPTS.py --inputfile '/path/to/file.txt' --time 300`

`python CalculPTS.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `inputfile` | Path to the input file (required) |
| `time`  | Time to make the calculations (required) |
#### Notes
The input file name should not contain a `.` except before the extension (*e.g.,* `input_file_name.sqlite`). The calculations only take into account the active infection. To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

## Calculate running times
This script calculates and summarizes the running times per runs and per replicates. It uses the `summary` and `meta` tables from the varmodel3 output database(s) (in [SQLite3 format](https://www.sqlite.org/fileformat.html)). The results will be stored into two output files: one for the running time per run and another one for the running times per replicate.
#### Example command
`python CalculRunTime.py --directory '/path/to/directory' --runs 20 --replicates 10`

`python CalculRunTime.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `directory` | Path to the directory containing the input file(s) (required) |
| `runs`  | Number of runs (required) |
| `replicates`  | Number of replicates per run (required) |
#### Notes
The input file(s) should be in a directory with the following structure: `directory/c<run>/r<replicate>` (*i.e.,* similar as the directory that is created when a parameter sweep is done). The two output files will have similar content if there is only one replicate per run. To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [statistics](https://docs.python.org/3/library/statistics.html), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

## Calculate prevalence
The script `CalculPreval.py` calculates the prevalence, the number of participants, and the number of hosts at a specific time.
It uses the "sampled_host" table from the varmodel3 output database (in [SQLite3 format](https://www.sqlite.org/fileformat.html)).
#### Example command
`python CalculPreval.py --inputfile '/path/to/file.txt' --time 300`

`python CalculPreval.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `inputfile` | Path to the input file (required) |
| `time`  | Time to make the calculations (days) (required) |

#### Notes
The input file name should not contain a `.` except before the extension (*e.g.,* `input_file_name.sqlite`).
The number of hosts and prevalence calculations only take into account the hosts with at least one active infection.
To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

## Calculate SNP call proportions
This script calculates the proportions of single nucleotide polymorphism (SNP) calls per host and per locus from a file in [THE REAL McCOIL categorical method format](https://github.com/EPPIcenter/THEREALMcCOIL/tree/master/categorical_method) ([Chang *et al.* 2017](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005348)). The proportion of heterozygote (`prop_hetero`; also known as double allele calls (DACs) or mixed allele calls (MACs)), homozygote major (`prop_homo_maj`), homozygote minor (`prop_homo_min`), and missing (`prop_miss`) calls are calculated. The results will be stored into two output files: 1) `file_name_PropCalls_Ind.txt` for the proportions of SNP calls per host; and 2) `file_name_PropCalls_Loc.txt` for the proportions of SNP calls per locus.
#### Example command
`python CalculPropCalls.py --inputfile '/path/to/file.txt`

`python CalculPropCalls.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `inputfile` | Path to the input file (required) |
#### Notes
The input file name should not contain a `.` except before the extension (*e.g.,* `input_file_name.txt`). To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

## Calculate targets
This script calculates the targets for the model fitting. These targets correspond to five diversity metrics: the prevalence, multiplicty of infection (MOI), pairwise type sharing (PTS), and number of strains and genes. It uses the 'sampled_host', 'sampled_infections', and 'sampled_infection_genes' tables from the varmodel3 output database (in [SQLite3 format](https://www.sqlite.org/fileformat.html)). When the `varmodel3` simulations involved very high diversity (*e.g.,* `n_genes_initial > 20000)`, the PTS is only calculated for a subsample of 1,000 strains to avoid potential memory issues. To account for measurement error, calculations can be done with a measurement model using `--measurement`. If so, it is required to provide an additional file describing the distribution of the number of gene per monoclonal infection ([Labbe *et al.* 2022](https://www.biorxiv.org/content/10.1101/2022.06.27.497801v1)). The later should contains two columns: 1) the number of genes per monoclonal infection, and 2) the weight reflecting the gene counts density function). If `--measurement` is used, it is also required to provide the proportion of host to keep for the calculations (e.g. proportion of infections that can be detected by microscopy).

#### Example command
`python Targets.py --inputfile '/path/to/file.txt' --time 300 --measurement --distribution '/path/to/distribution.txt' --prop 0.57`

`python Targets.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `inputfile` | Path to the input file (required) |
| `time` | Time to make the calculate the targets (required) |
| `measurement` | Wether calculations take into account a measurement model (optional) |
| `distribution` | Path to file describing the distribution of the number gene per monoclonal infection (optional) |
| `prop` | Proportion of host to keep in the calculations (optional) |
#### Notes
The input file name should not contain a `.` except before the extension (*e.g.,* `input_file_name.txt`). To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [numpy](https://numpy.org/), [random](https://docs.python.org/3/library/random.html), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

## Compare diversity metrics
This script compares the varmodel3 output databases (in [SQLite3 format](https://www.sqlite.org/fileformat.html)) located in two distinct directories. It could be used to evaluate new versions of the model, *e.g.,* outputs before vs. after a new implementation. Four major diversity metrics are compared: the average prevalence, pairwise type sharing (PTS), and number of strains and genes per replicate. For each diversity metric, the script plots its distributions and perform a [two-sample Kolmogorov-Smirnov test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.kstest.html). It uses the `sampled_hosts`, `gene_strain_counts`, `sampled_infections`, and `sampled_infection_genes` tables from the output database. 
#### Example command
`python DivMetComp.py --directory1 '/path/to/directory1' --directory2 '/path/to/directory2' --replicates 10`

`python DivMetComp.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `directory1` | Path to the directory containing the 1st set of varmodel3 output files (required) |
| `directory2` | Path to the directory containing the 2nd set of varmodel3 output files (required) |
| `replicates`  | Number of replicates per run (required) |
#### Notes
A minimum of 2 replicates is required, but we recommend using at least 10 replicates. To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [numpy](https://numpy.org/), [scipy](https://scipy.org/), [matplotlib](https://matplotlib.org/), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

## Convert SNPs in THE REAL McCOIL format
Based on the SNP allele frequencies, this script converts the SNP data into THE REAL McCOIL input format. See the `README.docx` file of [THE REAL McCOIL github](https://github.com/EPPIcenter/THEREALMcCOIL) for more details. SNP calling information is stored in a matrix, where each element Sij represents SNP information at locus j of individual i, and can be 0 [homozygous minor allele], 0.5 [heterozygous], 1 [homozygous major allele] or -1 [missing data]. The script also exports the SNP minor allele frequencies (MAF) in another output file. It uses the `sampled_infections` and `sampled_infection_snps` tables from the varmodel3 output database (in SQLite3 format).
#### Example command
`python TheRealMcCoilFormat.py --inputfile '/path/to/file.txt' --time 300 --minfreq 0.1`

`python TheRealMcCoilFormat.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `inputfile` | Path to the input file (required) |
| `time` | Time to make the calculations (required) |
| `minfreq`  | Minor allele frequency (MAF) for a SNP to be considered (required) |
#### Notes
The input file name should not contain a `.` except before the extension (*e.g.,* `input_file_name.txt`). To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [argparse](https://docs.python.org/3/library/argparse.html), and [sys](https://docs.python.org/3/library/sys.html).

## SNP measurement model
This script sample a number of missing loci and randomly replace them with missing data in an input file in [THE REAL McCOIL format](https://github.com/EPPIcenter/THEREALMcCOIL). It uses a distribution of the proportion of missing SNP loci per host based on Taq-Man assay empirical data. The distribution file should only contains two columns, i.e. the proportion and the weigths of missing SNP loci per host.
#### Example command
`python SnpErrMod.py --inputfile '/path/to/file.txt' --missdist '/path/to/file.txt'`

`python SnpErrMod.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `inputfile` | Path to the input file in THE REAL McCOIL format (required) |
| `missdist` | Path to the file providing the distribution of the proportion of missing loci per host (required) |
#### Notes
The input file name should not contain a `.` except before the extension (*e.g.,* `input_file_name.txt`). To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [argparse](https://docs.python.org/3/library/argparse.html), [random](https://docs.python.org/3/library/random.html), and [numpy](https://numpy.org/).
