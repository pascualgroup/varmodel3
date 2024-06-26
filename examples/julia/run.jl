#!/usr/bin/env julia

"""
The purpose of this file is to illustrate how to do a single ad-hoc run in Julia
without an external parameters file. This approach is convenient for testing.

Instead of using the standard `run.jl`, which loads parameters from JSON, you
generate parameters inside Julia and run the model code directly with them.

To use this script, do the following:

1. Copy this file into a new directory
2. Modify the relative paths to `preamble.jl` and `model.jl` below
3. Modify parameter values
4. Run the script: `julia run.jl`, or directly as a shell script, `./run.jl`
"""

# Load packages and definition of Params struct
include("../../preamble.jl")

# Define the parameters variable P, which exposes parameters to be used as
# compile-time constants when loading the model code below.
const P = let
    t_year = 360
    daily_biting_rate_multiplier = readdlm("../mosquito_population.txt", Float64)[:,1]
    #snp_ld_matrix = readdlm("../pairwise_ld_coefficient_24snps.txt", Float64)

    t_end_years = 100
    t_end = t_end_years * t_year

    # Uncomment this, and argument to Params() below, to enable an intervention
    # for some subset of years.
#     biting_rate_multiplier_by_year = repeat([1.0], t_end_years)
#     biting_rate_multiplier_by_year[61:62] .= 0.5

    t_burnin_years = 80
    t_burnin = t_burnin_years * t_year

    add_params(Params(), (
        upper_bound_recomputation_period = 30,

        output_db_filename = "output.sqlite",

        summary_period = 30,
        gene_strain_count_period = t_year,

        host_sampling_period = [30],
        host_sample_size = 100,

        verification_period = t_year,

        sample_infection_duration_every = 1000,

        rng_seed = nothing,

        whole_gene_immune = false,

        t_year = t_year,
        t_end = t_end,

        t_burnin = t_burnin,

        n_hosts = 10000,
        n_initial_infections = 20,

        n_genes_initial = 9600,
        n_genes_per_strain = 60,

        n_loci = 2,

        n_alleles_per_locus_initial = 960, 

        transmissibility = 1,
        coinfection_reduces_transmission = true,

        # ectopic_recombination_rate = 1.8e-7,
        ectopic_recombination_rate = [4.242641e-4, 4.242641e-4],
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

        switching_rate = [1.0/6.0, 1.0/6.0],
        # switching_rate = 1.0/6.0,

        mean_host_lifetime = 23.38837487739662 * t_year,

        background_clearance_rate = 0.0,

        immigration_rate_fraction = 0.0026,

        n_infections_liver_max = 20,
        n_infections_active_max = 20,

        biting_rate = 0.00002 * daily_biting_rate_multiplier, # 0.0005

#         biting_rate_multiplier_by_year = biting_rate_multiplier_by_year,

        migrants_match_local_prevalence = true,
        migration_rate_update_period = 30,


        
        # parameters for var groups implementation
        var_groups_functionality = [1, 1],
        var_groups_ratio = [0.25, 0.75],
        var_groups_fix_ratio = true,
        var_groups_do_not_share_alleles = true,
        var_groups_high_functionality_express_earlier = true,
        gene_group_id_association_recomputation_period = 30,
    ))
end

# Load and run the model code, which can now be compiled with reference to
# parameters constant P.
include("../../src/model.jl")
run()
