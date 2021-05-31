using JSON3
using Parameters
using Base.Filesystem
using StructTypes

"""
    Enum definition for the `model` parameter, described below.
"""
@enum Model begin
    VAR_ONLY = 1
    VAR_WITH_PARASITEMIA = 2
end

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
        (It it is possible that some of those could vary while preserving reproducibility.)
    """
    rng_seed::Union{Int, Missing} = missing
    
    """
        Which model to use.
        
        `VAR_ONLY` specifies the old model that does not include parasitemia, 
        which is close to the one implemented in the `varmodel2` C++ project.
    
        `VAR_WITH_PARASITEMIA` specifies the a model that combines the var-expression model
        with a parasitemia model based on OpenMalaria.
    
        The models are different enough that they are implemented separately, in the
        `var_only` and `var_with_parasitemia` directories, respectively.
    
        Both of these models have a number of switches to control behavior, as
        described below.
    
        In order to use parameter values as JIT-compile-time constants, rather than
        dynamically dispatching the appropriate version of the model, at compile
        time a different set of source files is read (see `model.jl`).
    """
    model::Union{Model, Missing} = missing
    
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
    upper_bound_recomputation_period::Union{Int, Missing} = missing
    
    """
        Filename for parasitemia patient data.
        
        Relative to working directory. Therefore, absolute path is recommended
        for sweeps.
    """
    patient_data_filename::Union{String, Missing} = missing
    
    """
        Filename for output database.
        
        Relative to working directory.
    """
    output_db_filename::Union{String, Missing} = missing
    
    """
        How often to sample host output.
        
        To turn off host sampling, set this to `missing` (`null` in JSON).
    """
    host_sampling_period::Union{Int, Missing} = missing
    
    """
        Number of hosts to sample at each sampling period.
    """
    host_sample_size::Union{Int, Missing} = missing
    
    
    """
        How often to write summary output.
    """
    summary_period::Union{Int, Missing} = missing
    
    
    """
        How often to output the number of circulating genes and strains.
        
        Because counts are not kept update dynamically at present, this may be
        the bottleneck in the simulation.
    """
    gene_strain_count_period::Union{Int, Missing} = missing
    
    """
        How often to verify consistency of simulation state.
        
        Useful primarily during development to find new bugs.
        If set to `missing` (`null` in JSON), no verification will be performed.
    """
    verification_period::Union{Int, Missing} = missing
    
    """
        Number of time units in a year.
        
        Currently used only to constrain size of `biting_rate_multiplier`.
    """
    t_year::Union{Int, Missing} = missing
    
    """
        Simulation end time.
    """
    t_end::Union{Int, Missing} = missing
    
    """
        Number of hosts.
        
        Currently constant through a simulation; deaths are coupled to births.
    """
    n_hosts::Union{Int, Missing} = missing
    
    """
        Number of initial infections.
        
        `n_initial_infections` hosts are each given 1 initial infection at t = 0.
    """
    n_initial_infections::Union{Int, Missing} = missing
    
    """
        Biting rate for each day of the year.
        
        Dimensions: (t_year,)
    """
    biting_rate::Union{Missing, Array{Float64}} = missing
    
    """
        Number of genes in strain.
        
        During an infection, each gene is expressed once, unless the host
        is already immune.
    """
    n_genes_per_strain::Union{Int, Missing} = missing
    
    """
        Number of genes in the initial gene pool.
        
        The gene pool is used to seed initial infections as well as infections
        due to immigration.
    """
    n_genes_initial::Union{Int, Missing} = missing
    
    """
        Number of epitope loci in each gene.
        
        A gene is comprised of a sequence of epitope alleles, one at each locus.
    """
    n_loci::Union{Int, Missing} = missing
    
    """
        Initial number of alleles for each epitope locus.
        
        The number of alleles changes over time due to mutation events.
    """
    n_alleles_per_locus_initial::Union{Int, Missing} = missing
    
    """
        Baseline transmissibility of infections.
        
        This is the probability that an infection in the source host is transmitted
        to the destination host during a biting event, if the source host
        has exactly one infection.
    """
    transmissibility::Union{Float64, Missing} = missing
    
    """
        Whether or not transmissibility is reduced with coinfection.
        
        If `true`, transmissibility is inversely proportional to the number
        of infections in the source host.
    """
    coinfection_reduces_transmission::Union{Bool, Missing} = missing
    
    """
        Ectopic recombination rate parameter.
        
        The number of recombinations, per active infection, per unit time
        is equal to:
        
        ```
        ectopic_recombination_rate * n_genes_per_strain * (n_genes_per_strain - 1) / 2
        ```
    """
    ectopic_recombination_rate::Union{Float64, Missing} = missing
    
    """
        Maximum immunity level.
        
        See description in the `immunity` field of struct `State`.
    """
    immunity_level_max::Union{Int8, Missing} = missing
    
    """
        Rate at which immunity is lost, per host, per gene.
        
        See description in the `immunity` field of struct `State`.
    """
    immunity_loss_rate::Union{Float64, Missing} = missing
    
    """
        Rate of mutation, per active infection.
        
        Mutation events randomly change a single epitope in a randomly chosen
        gene to a newly created allele.
    """
    mutation_rate::Union{Float64, Missing} = missing
    
    """
        Duration of the liver stage.
        
        Infections become active after `t_liver_stage` units of time.
    """
    t_liver_stage::Union{Float64, Missing} = missing
    
    """
        Switching rate for genes the host is not immune to.
        
        If the host is immune to the gene, then switching is instantaneous.
    """
    switching_rate::Union{Float64, Missing} = missing
    
    """
        Mean of exponential distribution used to draw host lifetime.
        
        This distribution is truncated at `max_host_lifetime`, so the effective
        mean is lower than this number.
    """
    mean_host_lifetime::Union{Float32, Missing} = missing
    
    """
        Maximum host lifetime.
        
        The exponential distribution used to draw host lifetime is truncated
        at this value.
    """
    max_host_lifetime::Union{Float32, Missing} = missing
    
    """
        Immigration rate, as a fraction of the non-immigration biting rate.
    """
    immigration_rate_fraction::Union{Float64, Missing} = missing
    
    """
        Maximum number of simultaneous infections in the liver stage.
        
        If a host is chosen for biting, the number of infections transmitted
        will be bounded by this number minus the number of current
        infections in the liver.
    """
    n_infections_liver_max::Union{Int, Missing} = missing
    
    """
        Maximum number of simultaneous active infections.
        
        When an infection is due to become active (from the liver stage), it is
        dropped if there are already `n_infections_active_max` infections.
    """
    n_infections_active_max::Union{Int, Missing} = missing
end

"""
    Instructs the JSON3 package to treat the Params type as simple data for the
    purposes of JSON serialization.
"""
StructTypes.StructType(::Type{Params}) = StructTypes.Struct()

"""
    
"""
function validate(p::Params)
    @assert p.model !== missing
    
    @assert p.upper_bound_recomputation_period !== missing
    @assert p.upper_bound_recomputation_period > 0
    
    @assert p.output_db_filename !== missing
    @assert p.output_db_filename != ""
    output_dir = dirname(p.output_db_filename)
    @assert output_dir == "" || isdir(output_dir)
    
    @assert p.host_sampling_period !== missing
    @assert p.host_sampling_period > 0
    @assert p.host_sample_size !== missing
    @assert p.host_sample_size >= 0
    
    @assert p.gene_strain_count_period !== missing
    @assert p.gene_strain_count_period > 0
    
    if p.verification_period !== missing
        @assert p.verification_period > 0
    end
    
    @assert p.t_year !== missing
    @assert p.t_year > 0
    
    @assert p.t_end !== missing
    @assert p.t_end >= 0
    
    @assert p.n_hosts !== missing
    @assert p.n_hosts > 0
    
    @assert p.n_initial_infections !== missing
    @assert p.n_initial_infections >= 0
    
    @assert p.biting_rate !== missing
    @assert length(p.biting_rate) == p.t_year
    
    @assert p.n_genes_initial !== missing
    @assert p.n_genes_initial > 0
    
    @assert p.n_genes_per_strain !== missing
    @assert p.n_genes_per_strain > 0
    
    @assert p.n_loci !== missing
    @assert p.n_loci > 0
    @assert p.n_loci == 2 # If different, need to revisit ectopic recombination model
    
    @assert p.n_alleles_per_locus_initial !== missing
    @assert p.n_alleles_per_locus_initial > 0
    
    @assert p.transmissibility !== missing
    @assert 0.0 <= p.transmissibility  <= 1.0
    
    @assert p.coinfection_reduces_transmission !== missing
    
    @assert p.ectopic_recombination_rate !== missing
    @assert p.ectopic_recombination_rate >= 0.0
    
    @assert p.immunity_loss_rate !== missing
    @assert p.immunity_loss_rate >= 0.0
    
    @assert p.mutation_rate !== missing
    @assert p.mutation_rate >= 0.0
    
    @assert p.t_liver_stage !== missing
    @assert p.t_liver_stage >= 0.0
    
    @assert p.switching_rate !== missing
    @assert p.switching_rate >= 0.0
    
    @assert p.mean_host_lifetime !== missing
    @assert p.mean_host_lifetime >= 0.0
    
    @assert p.max_host_lifetime !== missing
    @assert p.max_host_lifetime >= p.mean_host_lifetime
    
    @assert p.immigration_rate_fraction !== missing
    @assert p.immigration_rate_fraction >= 0.0
    
    @assert p.n_infections_liver_max !== missing
    @assert p.n_infections_liver_max >= 0
    
    @assert p.n_infections_active_max !== missing
    @assert p.n_infections_active_max >= 0
end
