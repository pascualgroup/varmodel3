# varmodel3

Ed Baskerville<br>
Frédéric Labbé

Var gene evolution model(s), implemented in Julia.

Based on previous C++ implementations (`varmodel` and `varmodel2`) by Ed Baskerville and Qixin He.

Model details are described inline in comments in the code; see *Code Organization* below to get oriented.

## Quickstart

### Setup

Before doing anything, run

```julia
./install-packages.jl
```

to install packages required by this code.

### Single runs

To do a single ad-hoc run of the model directly in Julia, copy the script `examples/julia/run.jl` into a new experiment directory, and modify and run as described in the comment string.

To perform a run with an existing parameters file in JSON format, copy the parameters file into a new experiment directory, and use the script `varmodel3/run.jl` as described in the comment string. To see how to generate parameters from JSON, see `examples/json/generate-json.jl`.

### Parameter sweeps

To do a parameter sweep, copy the `examples/sweep` directory, and modify/run as described in the comment string in `generate-sweep.jl`.

### Parameters
| Name | Description |
| :--: | ----------- | 
| `biting_rate` | Transmission rate for each day of the year |
| `coinfection_reduces_transmission`  | Whether or not transmissibility is reduced with coinfection |
| `distinct_initial_snp_allele_frequencies` | Whether the initial allele frequencies of the SNPs are distinct |
| `ectopic_recombination_generates_new_alleles` | Whether or not ectopic recombination generates new alleles |
| `ectopic_recombination_rate` | Ectopic recombination rate parameter |
| `gene_strain_count_period` | How often to output the number of circulating genes and strains |
| `host_sample_size` | Number of hosts to sample at each sampling period |
| `host_sampling_period` | How often to sample host output |

## Model Overview

TODO

## Code Organization

TODO
