#!/usr/bin/env julia

"""
The purpose of this file is to illustrate how to generate JSON from Julia for
running experiments.

This is useful for parameter sweeps, where you'll want to create an output
directory for each run, with a parameters file in each directory.

For a single ad-hoc run, there's no particular reason to generate JSON; see
`examples/julia` to see how to run the model directly.

To see this example expanded into a parameter sweep with multiple directories,
see `examples/sweep`.

To use this script, you'll want to copy this file to a different directory,
modify parameter settings, change the relative path to `preamble.jl` below,
and then do:

```
./generate-json.jl
<path-to>/varmodel3/run.jl parameters.json
```
"""

import JSON3
using DelimitedFiles

include("../../preamble.jl")

function main()
    params = init_params()
    validate(params)
    open("parameters.json", "w") do f
        JSON3.pretty(f, JSON3.write(params))
        println(f)
    end
end

function init_params()
    t_year = 360
    daily_biting_rate_multiplier = readdlm("../mosquito_population.txt", Float64)[:,1]

    Params(
        upper_bound_recomputation_period = 30,

        output_db_filename = "output.sqlite",

        summary_period = 30,
        gene_strain_count_period = 360,

        host_sampling_period = 30,
        host_sample_size = 100,

        verification_period = 360,

        rng_seed = missing,

        t_year = t_year,
        t_end = (111) * t_year,

        t_burnin = 61 * t_year,

        n_hosts = 10000,
        n_initial_infections = 20,

        n_genes_initial = 9600,
        n_genes_per_strain = 60,

        n_loci = 2,

        n_alleles_per_locus_initial = 960,

        transmissibility = 0.5,
        coinfection_reduces_transmission = true,

        ectopic_recombination_rate = 1.8e-7,

        immunity_level_max = 100,
        immunity_loss_rate = 0.001,

        mutation_rate = 1.42e-8,

        t_liver_stage = 14.0,

        switching_rate = 1.0/6.0,

        mean_host_lifetime = 30 * t_year,
        max_host_lifetime = 80 * t_year,

        immigration_rate_fraction = 0.0026,

        n_infections_liver_max = 10,
        n_infections_active_max = 10,

        biting_rate = 0.0005 * daily_biting_rate_multiplier,
    )
end

main()
