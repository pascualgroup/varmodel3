# Post analysis

This is a collection of scripts to analyze the sqlite output.
Below are notes about some useful tools. 
Not everything is documented yet, but most scripts have some helpful information if you type `python script.py -h`.

## Contents

* [Calculate running times](#Calculate-running-times)
* [Calculate prevalence](#Calculate-prevalence)

## Calculate running times
This script calculates and summarizes the running times per runs and per replicates. It uses the "summary" and "meta" tables from the varmodel3 output database(s) (in [SQLite3 format](https://www.sqlite.org/fileformat.html)). The results will be stored into two output files: one for the running time per run and another one for the running times per replicate.
#### Example command
`python CalculRunTime.py --directory '/path/to/directory' --runs 20 --replicates 10`

`python CalculRunTime.py -h` Will print a full list of command arguments.
#### Command arguments
| Name | Description |
| :--: | :---------: | 
| `directory` | Path to the directory containing the input file(s) |
| `runs`  | Number of runs |
| `replicates`  | Number of replicates per run |
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
| `inputfile` | Path to the input file |
| `time`  | Time to make the calculations (days) |

#### Notes
The input file name should not contains a "." except before the extension (e.g. "input_file_name.sqlite").
The number of hosts and prevalence calculations only take into account the hosts with at least one active infection.
To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: [os](https://docs.python.org/3/library/os.html), [sqlite3](https://docs.python.org/3/library/sqlite3.html), [pandas](https://pandas.pydata.org/), and [argparse](https://docs.python.org/3/library/argparse.html).
