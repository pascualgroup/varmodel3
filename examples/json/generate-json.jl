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

include("../../preamble.jl")

function main()
    params = init_params()
    validate(params)
    params_json = pretty_json(params)
    open("parameters.json", "w") do f
        println(f, params_json)
    end
end

function pretty_json(params)
    d = Dict(fn => getfield(params, fn) for fn in fieldnames(typeof(params)))
    io = IOBuffer()
    JSON.print(io, d, 2)
    String(take!(io))
end

function init_params()
    t_year = 360
    daily_biting_rate_multiplier = readdlm("../mosquito_population.txt", Float64)[:,1]
    #snp_ld_matrix = readdlm("../pairwise_ld_coefficient_24snps.txt", Float64)

    Params(
        upper_bound_recomputation_period = 30,

        output_db_filename = "output.sqlite",

        summary_period = 30,
        gene_strain_count_period = 360,

        host_sampling_period = 30,
        host_sample_size = 100,

        verification_period = 360,

        sample_infection_duration_every = 1000,

        rng_seed = nothing,

        whole_gene_immune = false,

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
        p_ectopic_recombination_is_conversion = 0.0,

        ectopic_recombination_generates_new_alleles = false,

#         ectopic_recombination_generates_new_alleles = true,
#         p_ectopic_recombination_generates_new_allele = 0.5,

        rho_recombination_tolerance = 0.8,
        mean_n_mutations_per_epitope = 5.0,

        immunity_level_max = 100,
        immunity_loss_rate = 0.001,

        mutation_rate = 1.42e-8,

        t_liver_stage = 14.0,

        switching_rate = 1.0/6.0,

        mean_host_lifetime = 30 * t_year,
        max_host_lifetime = 80 * t_year,

        background_clearance_rate = 0.0,

        immigration_rate_fraction = 0.0026,

        n_infections_liver_max = 10,
        n_infections_active_max = 10,

        biting_rate = 0.0005 * daily_biting_rate_multiplier,

        migrants_match_local_prevalence = true,
        migration_rate_update_period = 30,

        n_snps_per_strain = 24,

        distinct_initial_snp_allele_frequencies = false,
#         distinct_initial_snp_allele_frequencies = true,
#         initial_snp_allele_frequency = [0.1, 0.9],

        snp_linkage_disequilibrium = false,
#         snp_linkage_disequilibrium = true,
#         snp_pairwise_ld = snp_ld_matrix,
    )
end

main()
