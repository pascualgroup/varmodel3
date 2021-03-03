using JSON2
using Parameters
using Base.Filesystem
using StructTypes

@with_kw struct Params
    use_discrete_time_approximation::Union{Bool, Nothing} = nothing
    dt::Union{Int, Nothing} = nothing
    
    output_db_filename::Union{String, Nothing} = nothing
    
    output_hosts::Union{Bool, Nothing} = nothing
    output_strains::Union{Bool, Nothing} = nothing
    output_genes::Union{Bool, Nothing} = nothing
    
    host_sampling_period::Union{Int, Nothing} = nothing
    host_sample_size::Union{Int, Nothing} = nothing
    expected_equilibrium::Union{Int, Nothing} = nothing
    
    verification_on::Union{Bool, Nothing} = nothing
    verification_period::Union{Float64, Nothing} = nothing
    
    rng_seed::Union{UInt64, Nothing} = nothing
    
    t_year::Union{Int, Nothing} = nothing
    t_end::Union{Int, Nothing} = nothing
    t_burnin::Union{Int, Nothing} = nothing
    
    n_hosts::Union{Int, Nothing} = nothing
    n_initial_infections::Union{Int, Nothing} = nothing
    
    biting_rate_mean::Union{Nothing, Float64} = nothing
    daily_biting_rate_multiplier::Union{Nothing, Array{Float64}} = nothing
    
    n_genes_initial::Union{Int, Nothing} = nothing
    n_genes_per_strain::Union{Int, Nothing} = nothing
    
    n_loci::Union{Int, Nothing} = nothing
    
    n_alleles_per_locus_initial::Union{Int, Nothing} = nothing
    
    transmissibility::Union{Float64, Nothing} = nothing
    coinfection_reduces_transmission::Union{Bool, Nothing} = nothing
    
    ectopic_recombination_rate::Union{Float64, Nothing} = nothing
    
    immunity_level_max::Union{Int8, Nothing} = nothing
    immunity_loss_rate::Union{Float64, Nothing} = nothing
    
    mutation_rate::Union{Float64, Nothing} = nothing
    
    t_liver_stage::Union{Float64, Nothing} = nothing
    
    transition_rate_max::Union{Float64, Nothing} = nothing
    transition_rate_multiplier::Union{Array{Float64}, Nothing} = nothing
    
    mean_host_lifetime::Union{Float32, Nothing} = nothing
    max_host_lifetime::Union{Float32, Nothing} = nothing
    
    immigration_on::Union{Bool, Nothing} = nothing
    immigration_rate_fraction::Union{Float64, Nothing} = nothing
    
    n_infections_liver_max::Union{Int, Nothing} = nothing
    n_infections_active_max::Union{Int, Nothing} = nothing
end

StructTypes.StructType(::Type{Params}) = StructTypes.Struct()

function validate(p::Params)
    @assert p.use_discrete_time_approximation != nothing
    if p.use_discrete_time_approximation
        @assert p.dt != nothing
        @assert p.dt > 0
        @assert p.t_year % p.dt == 0
        @assert p.t_burnin % p.dt == 0
    end
    
    @assert p.output_db_filename != nothing
    @assert p.output_db_filename != ""
    output_dir = dirname(p.output_db_filename)
    @assert output_dir == "" || isdir(output_dir)
    
    @assert p.output_hosts != nothing
    @assert p.output_strains != nothing
    @assert p.output_genes != nothing
    
    if p.output_hosts
        @assert p.host_sampling_period != nothing
        @assert p.host_sampling_period > 0
        @assert p.host_sample_size >= 0
    end
    
    if p.expected_equilibrium != nothing
        @assert p.expected_equilibrium > 0
    end
    
    @assert p.verification_on != nothing
    if p.verification_on
        @assert p.verification_period != nothing
        @assert p.verification_period > 0
    end
    
    @assert p.t_year != nothing
    @assert p.t_year > 0
    
    @assert p.t_end != nothing
    @assert p.t_end >= 0
    @assert p.t_end % p.t_year == 0
    
    @assert p.t_burnin != nothing
    @assert p.t_burnin % p.t_year == 0
    
    @assert p.n_hosts != nothing
    @assert p.n_hosts >= 0
    
    @assert p.n_initial_infections != nothing
    @assert p.n_initial_infections >= 0
    
    @assert p.biting_rate_mean != nothing
    
    if p.daily_biting_rate_multiplier != nothing
        @assert length(p.daily_biting_rate_multiplier) == p.t_year
    end
    
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
    
    @assert p.transition_rate_max != nothing
    @assert p.transition_rate_max >= 0.0
    
    @assert p.transition_rate_multiplier != nothing
    @assert length(p.transition_rate_multiplier) == p.n_loci
    @assert all(0.0 .<= p.transition_rate_multiplier .<= 1.0)
    
    @assert p.mean_host_lifetime != nothing
    @assert p.mean_host_lifetime >= 0.0
    
    @assert p.max_host_lifetime != nothing
    @assert p.max_host_lifetime >= p.mean_host_lifetime
    
    @assert p.immigration_on != nothing
    
    @assert p.immigration_rate_fraction != nothing
    @assert p.immigration_rate_fraction >= 0.0
    
    @assert p.n_infections_liver_max != nothing
    @assert p.n_infections_liver_max >= 0
    
    @assert p.n_infections_active_max != nothing
    @assert p.n_infections_active_max >= 0
end
