# varmodel3

Ed Baskerville<br>
Frédéric Labbé

*Var* gene evolution model(s), implemented in Julia.
Based on previous C++ implementations ([varmodel](https://github.com/pascualgroup/varmodel) and [varmodel2](https://github.com/pascualgroup/varmodel2) by [Ed Baskerville](https://cobeylab.uchicago.edu/people/ed-baskerville/) and [Qixin He](https://www.qixinhe.net/?_ga=2.130167695.167565966.1656355847-504625709.1649177406). This code implements a model of malaria *var* gene evolution within an individual-based disease transmission model. Malaria strains are represented as unordered sets of *var* genes, which are in turn composed of abstract loci. A number of alleles can appear at each locus, and the allelic composition of a gene across loci governs immune dynamics in the host. Individual hosts are infected by strains, and infections can be transmitted between hosts. Each infection expresses a single *var* gene at a time, and the sequence of expressions is explicitly represented in the simulation. The simulation also includes immigration of new strains into the population, recombination during transmission and during an infection, and mutation. The simulation is modeled as a sequence of discrete events (state changes) that happen in continuous time. Model details are described inline in comments in the code; see *Code Organization* below to get oriented.

## Contents

* [Quickstart](#Quickstart)
* [Model overview](#Model-overview)
* [History of changes](#History-of-changes)
* [Parameters](#Parameters)
* [Code organization](#Code-organization)
* [Output database](#Output-database)
* [Citation](#Citation)

___
## Quickstart

### Setup

Before doing anything, run the following command line to install the packages required by this code:
```julia
./install-packages.jl
```

### Single runs

* A single ad-hoc run of the model, directly in Julia and without an external parameters file, is convenient for testing. Instead of using the standard `run.jl`, which loads parameters from JSON, you generate parameters inside Julia and run the model code directly with them. To do a single ad-hoc run of the model, do the following:
  1. Copy the script `examples/julia/run.jl` into a new directory,
  2. Modify the relative paths to `preamble.jl` and `model.jl`,
  3. Modify the parameter values,
  4. Run the script: `julia run.jl`, or directly as a shell script, `./run.jl`.

* To perform a run with an existing parameters file in JSON format, copy the parameters file into a new experiment directory, and use the script `varmodel3/run.jl` as described in the comment string. To see how to generate parameters from JSON, see `examples/json/generate-json.jl`. The purpose of this file is to illustrate how to generate JSON from Julia for running experiments. This is useful for parameter sweeps, where you will want to create an output directory for each run, with a parameters file in each directory.

### Parameter sweeps

To do a parameter sweep, copy the `examples/sweep` directory, and modify/run as described in the comment string in `generate-sweep.jl`. This script loops through parameter combinations, and replicates with different random seeds, and generates files necessary to perform runs on a local machine or on a SLURM cluster. It also divides runs into jobs suitable for execution on a single cluster node or local machine. The runs are specified as lines in the job's `runs.txt` file, and the job is specified in a `job.sbatch` file, which can be run directly as a shell script or submitted to a SLURM cluster. Each job uses the script `varmodel3/runmany.jl` to run a single-node, multi-core queue of runs, with one run running on each core at any time. This script also generates a script `submit_jobs.sh`, which submits every job to SLURM at once. Runs are divided into at most `N_JOBS_MAX` jobs that make use of at most `N_CORES_PER_JOB_MAX` for the cluster node's local queue. This allows you to work within limits set by your cluster administrator. If you have no limits, you should set `N_JOBS_MAX` to a very large number, and set `N_CORES_PER_JOB_MAX = 1`, so that the cluster can dynamically
balance runs across cluster nodes as the experiment runs. To modify configuration settings for SLURM jobs, edit the template string in the `generate_jobs()` function.

___
## Model overview

In this model, hosts carry infections of different strains of the malaria parasite *Plasmodium falciparum*. Each parasite genome consists of a specific combination (*i.e.,* repertoire) of `n_genes_per_strain` *var* genes. Strain identity is defined by this repertoire independent of order. Although unlikely, the same *var* gene may occur multiple times in a strain. Each *var* gene itself is represented as a linear combination of `n_loci` epitopes, *i.e.,* parts of the molecule that act as antigens and are targeted by the immune system. At the outset, each locus i has one of `n_alleles_per_locus_initial[i]` possible values, indexed from 0 to `n_alleles_per_locus_initial[i] - 1`. Mutation events create new alleles, so the number of distinct alleles at each locus can increase over time. At any time, hosts may be infected multiple times by the same or different strains. The *var* genes in a repertoire are expressed sequentially and the infection ends when the whole repertoire is depleted. The order of expression is randomized distinctly for each infection. The duration of the active period of a *var* gene, and thus of the infection, is determined by the number of unseen epitopes. When a *var* gene is deactivated, the host adds the deactivated *var* gene epitopes to its immunity memory. Specific immunity toward a given epitope experiences a loss rate from host immunity memory, and re-exposure is therefore required to maintain it. The local population is open to immigration from the regional pool.

___
## History of changes

This code is a new implementation in [Julia](https://julialang.org/) of the malaria *var* gene evolution model which is based on previous C++ implementations (*i.e.,* [varmodel](https://github.com/pascualgroup/varmodel) and [varmodel2](https://github.com/pascualgroup/varmodel2)). The main changes from the previous implementation are as follows:

* While the previous implementations of the stochastic agent-based model (ABM) were adapted from the next-reaction method which optimizes the Gillespie first-reaction method, this implementation uses a simpler [Gillespie algorithm](https://www.sciencedirect.com/science/article/pii/0021999176900413).
* Our model extension allows us to keep track of the neutral part of each migrant parasite genome assembled by sampling one of the two possible alleles at each of a defined number of neutral bi-allelic SNPs.
* While the extended model can generate homogeneous initial SNP allele frequencies by sampling the migrant alleles with an identical probability from the regional pool (*i.e.,* 0.5), it can also generate distinct initial SNP allele frequencies by sampling the migrant alleles from the regional pool with distinct probabilities that sum up to one (*e.g.,* 0.2 and 0.8) and are randomly picked from a defined range (*e.g.,* [0.1-0.9]).
* Moreover, to generate the neutral part of a recombinant parasite and mimic meiotic recombination, which happens within the mosquito during the sexual reproduction stage of the parasite, a random allele is sampled for each bi-allelic SNP.
* Finally, to allow for linkage disequilibrium (LD) across the neutral part of the genome, neutral bi-allelic SNPs can be non-randomly associated and co-segregate as defined in a matrix of LD coefficients indicating the probability that pairs of linked SNPs will co-segregate during the meiotic recombination.

___
## Parameters

The parameter names should match the variables defined in `src/parameters.jl`, and the values should match the appropriate [type](https://docs.julialang.org/en/v1/manual/types/).

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
| `n_snps_per_strain` | Number of biallelic neutral single nucleotide polymorphims (SNPs) in strain |
| `p_ectopic_recombination_is_conversion` | Probability that an ectopic recombination is a conversion |
| `rho_recombination_tolerance` | Recombination tolerance, rho, Drummond et al |
| `rng_seed` | Seed for random number generator |
| `sample_duration` | Sample an infection duration every `sample_duration` infection(s) |
| `snp_linkage_disequilibrium` | Whether the SNPs (or some SNPs) are in linkage disequilibrium (LD) |
| `snp_pairwise_ld` | Pairwise linkage disequilibrium (LD) matrix |
| `summary_period` | How often to write summary output |
| `switching_rate` | Switching rate for genes the host is not immune to |
| `transmissibility` | Baseline transmissibility of infections |
| `t_liver_stage` | Duration of the liver stage |
| `t_burnin` | Burn-in time |
| `t_end` | Simulation end time |
| `t_year` | Number of time units in a year |
| `upper_bound_recomputation_period` | How often to recompute upper bounds for rejection sampling |
| `use_immunity_by_allele` | Immunity model |
| `verification_period` | How often to verify consistency of simulation state |
| `whole_gene_immune` | Whether a host gains immunity towards a gene if the host has seen all the alleles |

___
## Code organization

TODO

___
## Output database

The output database is in [SQLite3 format](https://www.sqlite.org/fileformat.html), which can be easily accessed from [R](https://www.r-project.org/) using the [RSQLite library](https://cran.r-project.org/web/packages/RSQLite/index.html), from [Python](https://www.python.org/) using the built-in [sqlite3 library](https://docs.python.org/3/library/sqlite3.html), or from [Matlab](https://www.mathworks.com/products/matlab.html) using the [mksqlite package](http://mksqlite.sourceforge.net/). We also recommend using [VisiData](https://www.visidata.org/) or the graphical [SQLite browser](https://sqlitebrowser.org/), especially while testing.

#### Current database tables include:

| Name | Description |
| :--: | ----------- | 
| `gene_strain_counts` | Number of circulating strains and *var* genes in all sampled hosts at different sampling times |
| `initial_snp_allele_frequencies` | Initial allele frequencies of each SNP locus |
| `meta` | Information related to the run (*e.g.,* elapsed time) |
| `sampled_duration` | Information related to the infection durations in all sampled hosts (*e.g.,* infection time) |
| `sampled_host` | Information related to the sampled hosts at different sampling times (*e.g.,* birth and death times) |
| `sampled_immunity` | Immunity level of all sampled hosts at different sampling times |
| `sampled_infection_genes` | Information related to the *var* genes involved in the infections of the sampled hosts (*e.g.,* allele ID) |
| `sampled_infection_snps` | Information related to the SNP loci involved in the infections of the sampled hosts (*e.g.,* allele ID) |
| `sampled_infections` | Information related to the infections of the sampled hosts at different sampling times (*e.g.,* infection ID) |
| `summary` | Summary at different sampling times (*e.g.,* number of infections, bites, and infected hosts) |

___
## Citation

Please cite this when using the model:
[Labbé F, He Q, Zhan Q, Tiedje KE, Argyropoulos DC, Tan MH, et al. Neutral vs. non-neutral genetic footprints of *Plasmodium falciparum* multiclonal infections. bioRxiv; 2022. p. 2022.06.27.497801. doi:10.1101/2022.06.27.497801](https://www.biorxiv.org/content/10.1101/2022.06.27.497801v1).
