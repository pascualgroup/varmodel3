#!/usr/bin/env julia

"""
The purpose of this file is to illustrate how to do a single ad-hoc run in Julia.

Instead of using the provided `run.jl`, which loads parameters from JSON, you
generate parameters inside Julia and run the model code directly with them.

```
./generate-json.jl
../../run.jl parameters.json
```
"""

using DelimitedFiles

include("../../parameters.jl")
include("../../varmodel3.jl")

function main()
    params = init_params()
    run(params)
end

function init_params()
    t_year = 360
    daily_biting_rate_multiplier = readdlm("../mosquito_population.txt", Float64)[:,1]
    
    Params(
        use_discrete_time_approximation = true,
        dt = 1,

        output_db_filename = "output.sqlite",
        
        summary_period = 30,
        strain_count_period = 360,
        
        host_sampling_period = 30,
        host_sample_size = 100,
        
        verification_period = 30,

        rng_seed = nothing,

        t_year = t_year,
        t_end = 111 * t_year,

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

        switching_rate_max = 1.0/6.0,
        switching_rate_multiplier = [0.5, 1.0],

        mean_host_lifetime = 30 * t_year,
        max_host_lifetime = 80 * t_year,
        
        immigration_rate_fraction = 0.0026,

        n_infections_liver_max = 10,
        n_infections_active_max = 10,

        biting_rate = 0.0005 * daily_biting_rate_multiplier,
    )
end

main()
