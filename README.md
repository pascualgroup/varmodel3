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
| `immigration_rate_fraction` | Immigration rate, as a fraction of the non-immigration biting rate |
| `immunity_level_max` | Maximum immunity level |
| `immunity_loss_rate` | Rate at which immunity is lost, per host, per gene |
| `initial_snp_allele_frequency` | Range of the possible initial frequencies for one of the two SNP alleles |
| `max_host_lifetime` | Maximum host lifetime |
| `mean_host_lifetime` | Mean of exponential distribution used to draw host lifetime |
| `mean_n_mutations_per_epitope` | Mean number of mutations per epitope for similarity calculation |
| `migrants_match_local_prevalence` | Whether the immigration rate needs to time the local infection rate |
| `mutation_rate` | Rate of mutation, per active infection |
| `n_alleles_per_locus_initial` | Initial number of alleles for each epitope locus |
| `n_genes_initial` | Number of genes in the initial gene pool |
| `n_genes_per_strain` | Number of genes in strain |
| `n_hosts` | Number of hosts |
| `n_infections_active_max` | Maximum number of simultaneous active infections |
| `n_infections_liver_max` | Maximum number of simultaneous infections in the liver stage |
| `n_initial_infections` | Number of initial infections |
| `n_loci` | Number of epitope loci in each gene |
## Model Overview

TODO

## Code Organization

TODO
