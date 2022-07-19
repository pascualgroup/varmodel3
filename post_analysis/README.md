# Post analysis

This is a collection of scripts to analyze the sqlite output. Below are notes about some useful tools. Not everything is documented yet, but most scripts have some helpful information if you type `python script.py -h`.

## Contents

* [Calculate MOIvar](#Calculate-MOIvar)
* [Calculate running times](#Calculate-running-times)
* [Calculate prevalence](#Calculate-prevalence)

## Calculate MOIvar
This script calculates the multiplicity of infection (MOI) per hosts using the *var*coding approach. The *var*coding approach (also termed *var* genotyping or *var* fingerprinting), employs the highly polymorphic sequences encoding the immunogenic DBLα domain of PfEMP1 (*Plasmodium falciparum* erythrocyte membrane protein 1), the major surface antigen of the blood stage of infection ([Rask *et al.* 2010](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000933)). The multigene family known as *var* encodes variants of this surface antigen which can reach tens of thousands of variants in endemic populations ([Tonkin-Hill *et al.* 2021](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009269)). The extensive diversity of the *var* gene family together with the very low percentage of *var* genes shared between parasites facilitate measuring MOI by amplifying, pooling, sequencing, and counting the number of DBLα types in a host.

For each host, this script also exports from the simulations the number of unique *var* genes and the true MOI. To account for *var* gene potential sampling errors, an optional measurement model is also available, which sub-sampled the number of *var* genes per strains according to a weight reflecting the *var* gene counts density function. The script uses the `sampled_infections`, `sample_hosts`, and `sampled_infection_genes` tables from the varmodel3 output database (in [SQLite3 format](https://www.sqlite.org/fileformat.html)).

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
The MOI calculations only take into account the active infection(s). The input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite"). To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [random](https://docs.python.org/3/library/random.html), [numpy](https://numpy.org/), and [argparse](https://docs.python.org/3/library/argparse.html).

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
The input file(s) should be in a directory with the following structure: `directory/c<run>/r<replicate>` (*i.e.,* similar as the directory that is created when a parameter sweep is done). The two output files will have similar content if there is only one replicate per run. To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), [statistics](https://docs.python.org/3/library/statistics.html), and [argparse](https://docs.python.org/3/library/argparse.html).

## Calculate prevalence
The script `CalculPreval.py` calculates the prevalence, the number of participants, and the number of hosts at a specific time.
It uses the "sampled_host" table from the varmodel3 output database (in [SQLite3 format](https://www.sqlite.org/fileformat.html) format).
#### Example command
`python CalculPreval.py --inputfile '/path/to/file.txt' --time 300`

`python CalculPreval.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `inputfile` | Path to the input file (required) |
| `time`  | Time to make the calculations (days) (required) |

#### Notes
The input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
The number of hosts and prevalence calculations only take into account the hosts with at least one active infection.
To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), and [argparse](https://docs.python.org/3/library/argparse.html).
