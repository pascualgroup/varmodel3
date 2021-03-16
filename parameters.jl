using JSON2
using Parameters
using Base.Filesystem
using StructTypes

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
    rng_seed::Union{UInt64, Nothing} = nothing
    
    """
        Whether or not to simulate using a discrete-time approximation.
    """
    use_discrete_time_approximation::Union{Bool, Nothing} = nothing
    
    """
        Timestep of discrete-time approximation.
    """
    dt::Union{Int, Nothing} = nothing
    
    """
        Filename for output database.
        
        Relative to working directory.
    """
    output_db_filename::Union{String, Nothing} = nothing
    
    """
        How often to sample host output.
        
        To turn off host sampling, set this to `nothing` (`null` in JSON).
    """
    host_sampling_period::Union{Int, Nothing} = nothing
    
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
    strain_count_period::Union{Int, Nothing} = nothing
    
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
    ectopic_recombination_rate::Union{Float64, Nothing} = nothing
    
    """
        Maximum immunity level.
        
        See description in the `immunity` field of struct `State`.
    """
    immunity_level_max::Union{Int8, Nothing} = nothing
    
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
        Duration of the liver stage.
        
        Infections become active after `t_liver_stage` units of time.
    """
    t_liver_stage::Union{Float64, Nothing} = nothing
    
    """
        Switching rate for genes the host is not immune to.
        
        If the host is immune to the gene, then switching is instantaneous.
    """
    switching_rate::Union{Float64, Nothing} = nothing
    
    """
        Mean of exponential distribution used to draw host lifetime.
        
        This distribution is truncated at `max_host_lifetime`, so the effective
        mean is lower than this number.
    """
    mean_host_lifetime::Union{Float32, Nothing} = nothing
    
    """
        Maximum host lifetime.
        
        The exponential distribution used to draw host lifetime is truncated
        at this value.
    """
    max_host_lifetime::Union{Float32, Nothing} = nothing
    
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
end

"""
    Instructs the JSON3 package to treat the Params type as simple data for the
    purposes of JSON serialization.
"""
StructTypes.StructType(::Type{Params}) = StructTypes.Struct()

"""
    
"""
function validate(p::Params)
    @assert p.use_discrete_time_approximation != nothing
    if p.use_discrete_time_approximation
        @assert p.dt != nothing
        @assert p.dt == 1 # Code currently assumes this
        @assert p.dt > 0
        @assert p.t_year % p.dt == 0
    end
    
    @assert p.output_db_filename != nothing
    @assert p.output_db_filename != ""
    output_dir = dirname(p.output_db_filename)
    @assert output_dir == "" || isdir(output_dir)
    
    if p.host_sampling_period != nothing
        @assert p.host_sampling_period > 0
        @assert p.host_sample_size >= 0
    end
    
    if p.verification_period != nothing
        @assert p.verification_period > 0
    end
    
    @assert p.t_year != nothing
    @assert p.t_year > 0
    
    @assert p.t_end != nothing
    @assert p.t_end >= 0
    
    @assert p.n_hosts != nothing
    @assert p.n_hosts >= 0
    
    @assert p.n_initial_infections != nothing
    @assert p.n_initial_infections >= 0
    
    @assert p.biting_rate != nothing
    @assert length(p.biting_rate) == p.t_year
    
    @assert p.n_genes_initial != nothing
    @assert p.n_genes_initial > 0
    
    @assert p.n_genes_per_strain != nothing
    @assert p.n_genes_per_strain > 0
    
    @assert p.n_loci != nothing
    @assert p.n_loci > 0
    
    @assert p.n_alleles_per_locus_initial != nothing
    @assert p.n_alleles_per_locus_initial > 0
    
    @assert p.transmissibility != nothing
    @assert 0.0 <= p.transmissibility  <= 1.0
    
    @assert p.coinfection_reduces_transmission != nothing
    
    @assert p.ectopic_recombination_rate != nothing
    @assert p.ectopic_recombination_rate >= 0.0
    
    @assert p.immunity_level_max != nothing
    
    @assert p.immunity_loss_rate != nothing
    @assert p.immunity_loss_rate >= 0.0
    
    @assert p.mutation_rate != nothing
    @assert p.mutation_rate >= 0.0
    
    @assert p.t_liver_stage != nothing
    @assert p.t_liver_stage >= 0.0
    
    @assert p.switching_rate != nothing
    @assert p.switching_rate >= 0.0
    
    @assert p.mean_host_lifetime != nothing
    @assert p.mean_host_lifetime >= 0.0
    
    @assert p.max_host_lifetime != nothing
    @assert p.max_host_lifetime >= p.mean_host_lifetime
    
    @assert p.immigration_rate_fraction != nothing
    @assert p.immigration_rate_fraction >= 0.0
    
    @assert p.n_infections_liver_max != nothing
    @assert p.n_infections_liver_max >= 0
    
    @assert p.n_infections_active_max != nothing
    @assert p.n_infections_active_max >= 0
end
