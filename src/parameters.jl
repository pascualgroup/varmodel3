using Parameters
using Base.Filesystem

"""
Parameters for a simulation.

The `@with_kw` macro, provided by the `Parameters` package, generates a
keyword constructor for the class.
"""
@with_kw struct Params
    """
    Seed for random number generator.

    Runs should be reproducible with the same seed
    on the same version of Julia, on the same operating system,
    on the same processor architecture.
    (It is possible that some of those could vary while preserving reproducibility.)
    """
    rng_seed::Union{Int, Nothing} = nothing

    """
    How often to recompute upper bounds for rejection sampling.

    Currently, this is related only to immunity loss events.
    At any time, there is an upper bound `n` on the number of immunities
    possessed by any host in the system.

    For simplicity and memory savings, during immunity loss hosts are
    sampled uniformly randomly. If the host has `m` immune memories, the
    sampled host loses an immunity with probability `m / n`.

    Between recomputation periods, the upper bound will only grow.
    Every `upper_bound_recomputation_period`, the upper bound is reset to
    the exact value by scanning the number of immunities in each host.

    If the code is modified to allow for unlimited growth in the number of
    infections, this logic will apply to infections as well.
    """
    upper_bound_recomputation_period::Union{Int, Nothing} = nothing

    """
    Filename for output database.

    Relative to working directory.
    """
    output_db_filename::Union{String, Nothing} = nothing

    """
    How often to sample host output.

    To turn off host sampling, set this to `nothing` (`null` in JSON).
    """
    host_sampling_period::Union{Vector{Int}, Nothing} = nothing
    
    """
    Number of hosts to sample at each sampling period.
    """
    host_sample_size::Union{Int, Nothing} = nothing

    """
    How often to write summary output.
    """
    summary_period::Union{Int, Nothing} = nothing

    """
    How often to output the number of circulating genes and strains.

    Because counts are not kept update dynamically at present, this may be
    the bottleneck in the simulation.
    """
    gene_strain_count_period::Union{Int, Nothing} = nothing

    """
    How often to verify consistency of simulation state.

    Useful primarily during development to find new bugs.
    If set to `nothing` (`null` in JSON), no verification will be performed.
    """
    verification_period::Union{Int, Nothing} = nothing

    """
    Number of time units in a year.

    Currently used only to constrain size of `biting_rate_multiplier`.
    """
    t_year::Union{Int, Nothing} = nothing

    """
    Simulation end time.
    """
    t_end::Union{Int, Nothing} = nothing

    """
    Burn-in time.

    Begin writing output starting at `t == t_burnin`.
    """
    t_burnin::Union{Int, Nothing} = nothing

    """
    Number of hosts.

    Currently constant through a simulation; deaths are coupled to births.
    """
    n_hosts::Union{Int, Nothing} = nothing

    """
    Number of initial infections.

    `n_initial_infections` hosts are each given 1 initial infection at t = 0.
    """
    n_initial_infections::Union{Int, Nothing} = nothing

    """
    Biting rate for each day of the year.

    Dimensions: (t_year,)
    """
    biting_rate::Union{Nothing, Array{Float64}} = nothing

    """
    Biting rate multiplier for each year (optional).
    Used to apply yearlong interventions that reduce the biting rate.

    `biting_rate_multiplier_by_year[1]` corresponds to the time interval
    `[0, t_year)`, and so on.

    Dimensions: (t_end / t_year, )
    """
    biting_rate_multiplier_by_year::Union{Array{Float64}, Nothing} = nothing


    """
    Number of genes in strain.

    During an infection, each gene is expressed once, unless the host
    is already immune.
    """
    n_genes_per_strain::Union{Int, Nothing} = nothing

    """
    Number of genes in the initial gene pool.

    The gene pool is used to seed initial infections as well as infections
    due to immigration.
    """
    n_genes_initial::Union{Int, Nothing} = nothing

    """
    Number of epitope loci in each gene.

    A gene is comprised of a sequence of epitope alleles, one at each locus.
    """
    n_loci::Union{Int, Nothing} = nothing

    """
    Initial number of alleles for each epitope locus.

    The number of alleles changes over time due to mutation events.
    """
    n_alleles_per_locus_initial::Union{Int, Nothing} = nothing

    """
    Baseline transmissibility of infections.

    This is the probability that an infection in the source host is transmitted
    to the destination host during a biting event, if the source host
    has exactly one infection.
    """
    transmissibility::Union{Float64, Nothing} = nothing

    """
    Whether or not transmissibility is reduced with coinfection.

    If `true`, transmissibility is inversely proportional to the number
    of infections in the source host.
    """
    coinfection_reduces_transmission::Union{Bool, Nothing} = nothing

    """
    Ectopic recombination rate parameter.

    The number of recombinations, per active infection, per unit time
    is equal to:

    ```
    ectopic_recombination_rate * n_genes_per_strain * (n_genes_per_strain - 1) / 2
    ```
    """
    # ectopic_recombination_rate::Union{Float64, Nothing} = nothing
    ectopic_recombination_rate::Union{Vector{Float64}, Nothing} = nothing

    """
    Probability that an ectopic recombination is a conversion.
    """
    p_ectopic_recombination_is_conversion::Union{Float64, Nothing} = nothing

    """
    Whether or not ectopic recombination generates new alleles.

    If `true`, then the similarity calculation used to determine the probability
    that genes are functional will use real-valued breakpoints, and a new allele
    will be generated at the breakpoint with probability
    `p_ectopic_recombination_generates_new_allele`.
    """
    ectopic_recombination_generates_new_alleles::Union{Bool, Nothing} = nothing

    """
    If `ectopic_recombination_generates_new_alleles` is `true, then this is the
    probability that a new allele will be generated at the breakpoint.
    """
    p_ectopic_recombination_generates_new_allele::Union{Float64, Nothing} = nothing

    """
    Recombination tolerance, rho, Drummond et al.
    """
    rho_recombination_tolerance::Union{Float64, Nothing} = nothing

    """
    Mean number of mutations per epitope for similarity calculation.
    """
    mean_n_mutations_per_epitope::Union{Float64, Nothing} = nothing

    """
    Maximum immunity level.

    See description in the `immunity` field of struct `State`.
    """
    immunity_level_max::Union{Int16, Nothing} = nothing

    """
    Rate at which immunity is lost, per host, per gene.

    See description in the `immunity` field of struct `State`.
    """
    immunity_loss_rate::Union{Float64, Nothing} = nothing

    """
    Rate of mutation, per active infection.

    Mutation events randomly change a single epitope in a randomly chosen
    gene to a newly created allele.
    """
    mutation_rate::Union{Float64, Nothing} = nothing

    """
    Mean duration of the liver stage.

    Infections become active after `t_liver_stage` units of time on average.
    Actual progression is governed by an Erlang distribution with shape `liver_erlang_shape`.
    """
    t_liver_stage::Union{Float64, Nothing} = nothing

    """
    Shape parameter for Erlang distribution governing liver stage progress.
    
    Default shape of 49 gives a standard deviation of 2 days with a mean duration of 14 days.
    """
    liver_erlang_shape::Union{Int, Nothing} = 49

    """
    Switching rate for genes the host is not immune to.

    If the host is immune to the gene, then switching is instantaneous.
    """
    switching_rate::Union{Vector{Float64}, Nothing} = nothing
    # switching_rate::Union{Float64, Nothing} = nothing

    """
    Mean of exponential distribution used to draw host lifetime.

    This distribution is truncated at `max_host_lifetime`, so the effective
    mean is lower than this number.
    """
    mean_host_lifetime::Union{Float32, Nothing} = nothing

    """
        Maximum host lifetime.

        Not used, should be `nothing`.
    """
    max_host_lifetime::Union{Float32, Nothing} = nothing

    """
        Background clearance rate due to processes not explicitly modeled.

        Per active infection, per unit time (day).
    """
    background_clearance_rate::Union{Float64, Nothing} = nothing

    """
        Immigration rate, as a fraction of the non-immigration biting rate.
    """
    immigration_rate_fraction::Union{Float64, Nothing} = nothing

    """
        Maximum number of simultaneous infections in the liver stage.

        If a host is chosen for biting, the number of infections transmitted
        will be bounded by this number minus the number of current
        infections in the liver.
    """
    n_infections_liver_max::Union{Int, Nothing} = nothing

    """
        Maximum number of simultaneous active infections.

        When an infection is due to become active (from the liver stage), it is
        dropped if there are already `n_infections_active_max` infections.
    """
    n_infections_active_max::Union{Int, Nothing} = nothing

    """
        Parameter controlling how often to record infection durations.

        Every `sample_infection_duration_every` clearances, duration and other
        related information will be recorded about the infection.

        Q: would the statistical distribution be affected if the infections were
        chosen based on infection time rather than clearance time?
    """
    sample_infection_duration_every::Union{Int, Nothing} = nothing

    """
        Whether a host gains immunity towards a gene if the host has seen all
        the alleles. If `true`, then a host gains immunity towards a gene only
        if the host has seen all the alleles. If `false`, then a host gains a
        gradual increase in immunity as the host sees more alleles in a gene.
    """
    whole_gene_immune::Union{Bool, Nothing} = nothing

    """
        Whether the immigration rate needs to time the local infection rate.
        If `true`, then the immigration rate will be multiplied by the local
        infection rate (i.e. number of transmitted bites / number of total bites).
        If `false`, then the immigration rate is not impacted by the local infection
        rate.
    """
    migrants_match_local_prevalence::Union{Bool, Nothing} = nothing

    """
        How often to update migration rate based on local prevalence.
    """
    migration_rate_update_period::Union{Int, Nothing} = nothing

    """
        Number of biallelic neutral single nucleotide polymorphims (SNPs) in strain.
        These SNPs do not contribute to the infection or host immune memory.
        They are used to keep track of the neutral part of each parasite genome.
    """
    n_snps_per_strain::Union{Int, Nothing} = nothing

    """
        Whether the initial allele frequencies of the SNPs are distinct.
        If `true`, then the initial allele frequency are distinct in each SNP.
        For example, while one SNP could have its minor allele frequency (MAF)
        equals to 0.5, another SNP could have its MAF equals to 0.3.
        If `false`, then the initial allele frequencies are similar in all SNPs,
        i.e. their initial MAFs equal 0.5.
    """
    distinct_initial_snp_allele_frequencies::Union{Bool, Nothing} = nothing

    """
        If `distinct_initial_snp_allele_frequencies` is `true, then this is the
        range of the possible initial frequencies for one of the two SNP alleles,
        e.g. [0.1, 0.9].
    """
    initial_snp_allele_frequency::Union{Array{Float32}, Nothing} = nothing

    """
        Whether the SNPs (or some SNPs) are in linkage disequilibrium (LD).
        If `true`, then the linked SNPs will have similar initial SNP allele
        frequencies and they will tend to co-segregate during recombination.
        If `false`, then all the SNPs are considered unlinked and they will
        evolve independently.
    """
    snp_linkage_disequilibrium::Union{Bool, Nothing} = nothing

    """
        If `snp_linkage_disequilibrium` is `true`, then this is the pairwise
        linkage disequilibrium (LD) matrix. This matrix should provide the
        coefficient of LD between each pair of SNPs. If two loci are not
        coinherited at all (they are independent) then the value will be 0.0.
        If two loci are in total disequilibrium then value would be 1.0.
        The R script "Rscript_Create_LDPairwiseMatrix.R" could be used to
        create a pairwise LD matrix.
    """
#     snp_pairwise_ld::Union{Array{Float32, 2}, Nothing} = nothing
    snp_pairwise_ld::Union{Vector{Vector{Float64}}, Nothing} = nothing

    """
        below is the parameters added for implementing different var groups 
        in this Julia code; also modified switching_rate paremeter above from 
        a single value to a vector with elements corresponding to different
        groups.
    """
    var_groups_functionality::Union{Vector{Float32}, Nothing} = nothing
    var_groups_ratio::Union{Vector{Float32}, Nothing} = nothing
    var_groups_fix_ratio::Union{Bool, Nothing} = nothing
    var_groups_do_not_share_alleles::Union{Bool, Nothing} = nothing
    var_groups_high_functionality_express_earlier::Union{Bool, Nothing} = nothing
    gene_group_id_association_recomputation_period::Union{Int, Nothing} = nothing
    
    """
        below is the additional parameters. 
    """
    irs_start_year::Union{Int, Nothing} = nothing
    irs_duration::Union{Int, Nothing} = nothing
    biting_rate_factor::Union{Float64, Nothing} = nothing
    t_host_sampling_start::Union{Int, Nothing} = nothing 
    biting_rate_mean::Union{Float64, Nothing} = nothing
    # t_decimal_advance::Union{Float64, Nothing} = nothing

    """
        Whether to gather a sampling performance profile from this run.
    """
    profile_on::Bool = false

    """
        Filename for profiling data
    """
    profile_filename::String = "profile_data"

    """
        Delay between samples
    """
    profile_delay::Float64 = 0.1

    """
        Batch size for batched exponential draws.
        Exponential draws are a significant percent of simulation time;
        drawing them in large batches enables the use of vector (SIMD) instructions
        for a ~ 4x RNG speedup in isolated tests.
    """
    rng_batch_size::Int = 10000000
end

"""

"""
function validate(p::Params)
    @assert p.upper_bound_recomputation_period !== nothing
    @assert p.upper_bound_recomputation_period > 0

    @assert p.output_db_filename !== nothing
    @assert p.output_db_filename != ""
    output_dir = dirname(p.output_db_filename)
    @assert output_dir == "" || isdir(output_dir)

    @assert all(p.host_sampling_period.!==nothing)
    @assert all(p.host_sampling_period.>=0)

    @assert p.host_sample_size !== nothing
    @assert p.host_sample_size >= 0

    @assert p.gene_strain_count_period !== nothing
    @assert p.gene_strain_count_period > 0

    @assert p.sample_infection_duration_every !== nothing
    @assert p.sample_infection_duration_every >= 0

    if p.verification_period !== nothing
        @assert p.verification_period > 0
    end

    @assert p.t_year !== nothing
    @assert p.t_year > 0

    @assert p.t_end !== nothing
    @assert p.t_end >= 0

    @assert p.t_burnin === nothing || p.t_burnin >= 0

    @assert p.n_hosts !== nothing
    @assert p.n_hosts > 0

    @assert p.n_initial_infections !== nothing
    @assert p.n_initial_infections >= 0

    @assert p.biting_rate !== nothing
    @assert length(p.biting_rate) == p.t_year

    if p.biting_rate_multiplier_by_year !== nothing
        @assert length(p.biting_rate_multiplier_by_year) == Int(p.t_end / p.t_year)
        @assert all(0.0 .<= p.biting_rate_multiplier_by_year .<= 1.0)
    end

    @assert p.n_genes_initial !== nothing
    @assert p.n_genes_initial > 0

    @assert p.n_genes_per_strain !== nothing
    @assert p.n_genes_per_strain > 0

    @assert p.n_loci !== nothing
    @assert p.n_loci > 0
    # @assert p.n_loci == 2 # If different, need to revisit ectopic recombination model
    @assert p.n_loci >= 1

    @assert p.n_alleles_per_locus_initial !== nothing
    @assert p.n_alleles_per_locus_initial > 0

    @assert p.transmissibility !== nothing
    @assert 0.0 <= p.transmissibility <= 1.0

    @assert p.coinfection_reduces_transmission !== nothing

    @assert p.ectopic_recombination_rate !== nothing
    # @assert p.ectopic_recombination_rate >= 0.0
    @assert all(p.ectopic_recombination_rate .>= 0.0)

    @assert !isnothing(p.p_ectopic_recombination_is_conversion)
    @assert 0.0 <= p.p_ectopic_recombination_is_conversion <= 1.0

    @assert !isnothing(p.ectopic_recombination_generates_new_alleles)
    if p.ectopic_recombination_generates_new_alleles
        @assert !isnothing(p.p_ectopic_recombination_generates_new_allele)
        @assert 0.0 <= p.p_ectopic_recombination_generates_new_allele <= 1.0
    else
        @assert isnothing(p.p_ectopic_recombination_generates_new_allele)
    end

    @assert !isnothing(p.rho_recombination_tolerance)
    @assert 0.7 <= p.rho_recombination_tolerance <= 0.9

    @assert !isnothing(p.mean_n_mutations_per_epitope)
    @assert p.mean_n_mutations_per_epitope > 0.0

    @assert p.immunity_loss_rate !== nothing
    @assert p.immunity_loss_rate >= 0.0

    @assert p.mutation_rate !== nothing
    @assert p.mutation_rate >= 0.0

    @assert p.t_liver_stage !== nothing
    @assert p.t_liver_stage >= 0.0

    @assert p.liver_erlang_shape !== nothing
    @assert p.liver_erlang_shape >= 1
    @assert p.liver_erlang_shape <= 100

    @assert p.switching_rate !== nothing
    @assert all(p.switching_rate .>= 0.0)
    # @assert p.switching_rate >= 0.0

    @assert p.mean_host_lifetime !== nothing
    @assert p.mean_host_lifetime >= 0.0

    @assert p.max_host_lifetime === nothing

    @assert p.background_clearance_rate !== nothing
    @assert p.background_clearance_rate >= 0.0

    @assert p.immigration_rate_fraction !== nothing
    @assert p.immigration_rate_fraction >= 0.0

    if p.n_infections_liver_max !== nothing
        @assert p.n_infections_liver_max >= 0
    end

    if p.n_infections_active_max !== nothing
        @assert p.n_infections_active_max >= 0
    end

    @assert p.whole_gene_immune !== nothing

    @assert p.migrants_match_local_prevalence !== nothing
    if p.migrants_match_local_prevalence
        @assert p.migration_rate_update_period !== nothing
    end

    @assert p.n_snps_per_strain !== nothing
    @assert p.n_snps_per_strain >= 0
    if p.n_snps_per_strain > 0
        @assert p.distinct_initial_snp_allele_frequencies !== nothing
        if p.distinct_initial_snp_allele_frequencies
            @assert p.initial_snp_allele_frequency !== nothing
            @assert length(p.initial_snp_allele_frequency) == 2
            @assert p.initial_snp_allele_frequency[1] >= 0.0
            @assert p.initial_snp_allele_frequency[1] < 1.0
            @assert p.initial_snp_allele_frequency[2] > 0.0
            @assert p.initial_snp_allele_frequency[2] <= 1.0
            @assert p.initial_snp_allele_frequency[1] < p.initial_snp_allele_frequency[2]
        end
        @assert p.snp_linkage_disequilibrium !== nothing
        if p.snp_linkage_disequilibrium
            @assert p.n_snps_per_strain >= 2
            @assert length(p.snp_pairwise_ld) == p.n_snps_per_strain
            for i in 1:p.n_snps_per_strain
                @assert length(p.snp_pairwise_ld[i]) == p.n_snps_per_strain
            end
        end
    end

    """
    check the parameters for different var groups implementation.
    """
    @assert p.var_groups_functionality !== nothing
    @assert all(p.var_groups_functionality .>= 0) && all(p.var_groups_functionality .<= 1)
    @assert p.var_groups_ratio !== nothing
    @assert all(p.var_groups_ratio .>= 0) && all(p.var_groups_ratio .<= 1)
    @assert p.var_groups_fix_ratio !== nothing
    @assert p.var_groups_do_not_share_alleles !== nothing
    @assert p.var_groups_high_functionality_express_earlier !== nothing
    @assert all(round.(p.var_groups_ratio * p.n_genes_initial) .== p.var_groups_ratio * p.n_genes_initial) # check the number of genes in each group is an integer number
    @assert all(round.(p.var_groups_ratio * p.n_alleles_per_locus_initial) .== p.var_groups_ratio * p.n_alleles_per_locus_initial)
    @assert p.gene_group_id_association_recomputation_period !== nothing
    @assert p.gene_group_id_association_recomputation_period > 0
    
    """
    check additional params
    """
    @assert p.irs_start_year === nothing || p.irs_start_year >= 0 
    @assert p.irs_duration === nothing || p.irs_duration >= 0 
    @assert p.biting_rate_factor === nothing || p.biting_rate_factor >= 0.0
    @assert p.t_host_sampling_start === nothing || p.t_host_sampling_start >= 0
    @assert p.biting_rate_mean !== nothing && p.biting_rate_mean > 0.0
    # @assert p.t_decimal_advance !== nothing && p.t_decimal_advance > 0.0
end
