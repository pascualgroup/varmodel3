# Post analysis

This is a collection of scripts to analyze the sqlite output.
Below are notes about some useful tools. 
Not everything is documented yet, but most scripts have some helpful information if you type `python script.py -h`.

## Contents

* [Calculate prevalence](#Calculate-prevalence)

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
To run it, ensure that you are using Python v.3.7, and have installed the following dependencies: os, sqlite3, pandas, and argparse.
