#!/usr/bin/env julia

"""
The purpose of this file is to illustrate how to do a parameter sweep in Julia.

This script loops through parameter combinations, and replicates with different
random seeds, and generates files necessary to perform runs on a local machine
or on a SLURM cluster.

For each run, it creates a directory, `<combo_id>/<replicate_id>`, and adds
entries to a SQLite database of run information, to make it easy to identify
runs and collate output.

It also divides runs into jobs suitable for execution on a single cluster node
(or local machine). Each job is specified by a file

To use this file for a real experiment, you should copy it to your experiment
directory, modify the parameter combinations, and then run:

./generate-runs.jl

to generate the files for running the experiment.

...TODO
"""

using SQLite
using DelimitedFiles

include("../../parameters.jl")
include("../../varmodel3.jl")

const N_REPLICATES = 10

const N_JOBS = 100
const N_CORES_PER_JOB = 28

function main()
    
    base_params = init_base_params()
    
    combo_id = 1
    for biting_rate_mean in (0.0005, 0.0010, 0.0015)
        for transmissibility in (0.25, 0.5, 0.75)
            for replicate_id in 1:N_REPLICATES
                
            end
            group_id += 1
        end
    end
end

function init_base_params()
    t_year = 360
    daily_biting_rate_multiplier = readdlm("../mosquito_population.txt", Float64)[:,1]
    
    Params(
        use_discrete_time_approximation = true,
        dt = 1,

        output_db_filename = "output.sqlite",

        output_hosts = true,
        output_strains = false,
        output_genes = false,

        host_sampling_period = 30,
        host_sample_size = 100,
        expected_equilibrium = 10800,

        verification_on = true,
        verification_period = 30,

        rng_seed = nothing,

        t_year = t_year,
        t_end = 111 * t_year,
        t_burnin = 10 * t_year,

        n_hosts = 10000,
        n_initial_infections = 20,

        n_genes_initial = 9600,
        n_genes_per_strain = 60,

        n_loci = 2,

        n_alleles_per_locus_initial = 960,

        transmissibility = 0.5,
        coinfection_reduces_transmission = true,

        ectopic_recombination_rate = 1.8e-7,

        max_immunity_count = 100,
        immunity_loss_rate = 0.001,

        mutation_rate = 1.42e-8,

        t_liver_stage = 14.0,

        transition_rate_max = 1.0/6.0,
        transition_rate_multiplier = [0.5, 1.0],

        mean_host_lifetime = 30 * t_year,
        max_host_lifetime = 80 * t_year,

        immigration_on = true,
        immigration_rate_fraction = 0.0026,

        max_infection_count = 9,

        biting_rate_mean = 0.0005,
        daily_biting_rate_multiplier = daily_biting_rate_multiplier
    )
end

main()
