#!/usr/bin/env julia

"""
The purpose of this file is to illustrate how to generate JSON from Julia for
running experiments.

This is useful for parameter sweeps, where you'll want to create an output
directory for each run, with a parameters file in each directory.

For a single ad-hoc run, there's no particular reason to generate JSON; see
`examples/single-julia` to see how to run the model directly.

To see this example expanded into a parameter sweep with multiple directories,
see `examples/sweep`.

To actually perform the run, you need to (1) run this script and then
(2) run the provided `run.jl` file with the output:

```
./generate-json.jl
../../run.jl parameters.json
```
"""

import JSON3
using DelimitedFiles

include("../../parameters.jl")

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
        implementation = DISCRETE_TIME,
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

        transition_rate = 1/6.0,

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
