using DataStructures
using Parameters
using Random
using Distributions
using StatsBase
using StaticArrays
using Dates
using Test

const StrainId = UInt32
const AlleleId = UInt16
const Locus = UInt8
const ExpressionIndex = UInt8
const ImmunityLevel = Int8
const InfectionCount = Int8

@with_kw mutable struct State
    "Number of alleles for each locus (epitope)"
    n_alleles::Vector{AlleleId}
    
    """
        Gene pool used to generate initial infections and immigration events.
        Entry (i, j) is the allele ID for gene i, locus j.
        
        Currently fixed-size.
        Dimensions: (n_genes_initial, n_loci)
    """
    gene_pool::Array{AlleleId, 2}
    
    """
        ID for next strain to be created.
        
        Whenever a new strain is created via infection, immigration, mutation,
        or ectopic recombination, it is given the next strain ID and this
        counter is incremented.
    """
    next_strain_id::StrainId
    
    """
        Birth time for each host.
        
        Dimensions: (n_hosts,)
    """
    t_birth::Vector{Float32}
    
    """
        Death time for each host.
        
        Dimensions: (n_hosts,)
    """
    t_death::Vector{Float32}
    
    """
        Immunity level for each host to each epitope allele at each locus.
        Entry (i, j, k) is the immunity level for host i, allele j, locus k.
        
        When a host is exposed to an allele, its immunity level is incremented
        by 1, up to a maximum of immunity_level_max.
        
        During immunity loss events, immunity levels are decremented by 1,
        down to a minimum of 0.
        
        The array has an entry for every allele in the system for every host;
        the size of the array is an upper bound on the number of alleles at
        each locus. When the number of alleles exceeds the upper bound,
        the size of the array is increased by 25%.
        
        Dimensions: (n_hosts, n_alleles_upper_bound, n_loci)
    """
    immunity::Array{ImmunityLevel, 3}
    
    """
        Time of infection for infections in the liver stage.
        
        Fixed-size array for the sake of predictable memory usage.
        Unused slots are indicated by NaN.
        
        Dimensions: (n_hosts, n_infections_liver_max)
    """
    t_infection_liver::Array{Float32, 2}
    
    """
        Time of infection for active infections.
        
        Fixed-size array for the sake of predictable memory usage.
        Unused slots are indicated by NaN.
        
        Dimensions: (n_hosts, n_infections_active_max)
    """
    t_infection_active::Array{Float32, 2}
    
    """
        Strain IDs for infections in the liver stage.
        
        Whenever a strain is created or modified, it is given a new strain ID.
        Unused slots are indicated by 0.
        
        Dimensions: (n_hosts, n_infections_liver_max)
    """
    strain_id_liver::Array{StrainId, 2}
    
    """
        Strain IDs for active infections, transferred directly from the liver stage.
        
        Unused slots are indicated by 0.
        
        Dimensions: (n_hosts, n_infections_active_max)
    """
    strain_id_active::Array{StrainId, 2}
    
    """
        Sequence of genes for infections in the liver stage.
        
        Stored directly as a matrix of allele IDs.
        
        Dimensions: (n_hosts, n_infections_liver_max, n_genes_per_strain, n_loci,)
    """
    genes_liver::Array{AlleleId, 4}
    
    """
        Sequence of genes for infections in the liver stage.
        
        Transferred directly from genes_liver upon activation.
        Stored directly as a matrix of allele IDs.
        
        Dimensions: (n_hosts, n_infections_active_max, n_genes_per_strain, n_loci,)
    """
    genes_active::Array{AlleleId, 4}
    
    """
        Index of currently expressed gene for each active infection.
        
        Dimensions: (n_hosts, n_infections_active_max,)
    """
    expression_index::Array{ExpressionIndex, 2}
end

function State(p::Params)
    lifetime = Vector([draw_host_lifetime(p) for i in 1:p.n_hosts])
    t_birth = -rand(Float32, p.n_hosts) .* lifetime
    t_death = t_birth + lifetime
    
#     n_infections_liver = fill(0, p.n_hosts)
#     n_infections_active = fill(0, p.n_hosts)
    
    t_infection_liver = fill(NaN32, p.n_hosts, p.n_infections_liver_max)
    t_infection_active = fill(NaN32, p.n_hosts, p.n_infections_active_max)
    strain_id_liver = fill(0, p.n_hosts, p.n_infections_liver_max)
    strain_id_active = fill(0, p.n_hosts, p.n_infections_active_max)
    genes_liver = fill(AlleleId(0), p.n_hosts, p.n_infections_liver_max, p.n_genes_per_strain, p.n_loci)
    genes_active = fill(AlleleId(0), p.n_hosts, p.n_infections_active_max, p.n_genes_per_strain, p.n_loci)
    expression_index = fill(ExpressionIndex(0), p.n_hosts, p.n_infections_active_max)
    
    gene_pool = reshape(rand(
        AlleleId(1):AlleleId(p.n_alleles_per_locus_initial),
        p.n_genes_initial * p.n_loci
    ), (p.n_genes_initial, p.n_loci))
    
    infection_hosts = sample(1:p.n_hosts, p.n_initial_infections, replace = false)
#     n_infections[infection_hosts] .= 1
#     n_infections_liver[infection_hosts] .= 1
    t_infection_liver[infection_hosts] .= 0.0
    
    strain_id_liver[infection_hosts, 1] = 1:p.n_initial_infections
    next_strain_id = p.n_initial_infections + 1
    
    for i in 1:p.n_genes_per_strain
        genes_liver[infection_hosts, 1, i, :] = gene_pool[rand(1:p.n_genes_initial, p.n_initial_infections), :]
    end
    
    State(
        t_birth = t_birth,
        t_death = t_death,
        
#         n_infections = n_infections,
        t_infection_liver = t_infection_liver,
        t_infection_active = t_infection_active,
        
        strain_id_liver = strain_id_liver,
        strain_id_active = strain_id_active,
        
        next_strain_id = next_strain_id,
        
        genes_liver = genes_liver,
        genes_active = genes_active,
        
        expression_index = expression_index,
        
        n_alleles = fill(AlleleId(p.n_alleles_per_locus_initial), p.n_loci),
        immunity = fill(0, p.n_hosts, p.n_alleles_per_locus_initial, p.n_loci),
        
        gene_pool = gene_pool
    )
end

function verify(p::Params, s::State)
    @assert all(.!isnan.(s.t_birth))
    @assert all(.!isnan.(s.t_death))
    
    # Check infection times; use isnan to identify liver-stage and active infections
    liver_indices = findall(.!isnan.(s.t_infection_liver))
    liver_indices_null = findall(isnan.(s.t_infection_liver))
    active_indices = findall(.!isnan.(s.t_infection_active))
    active_indices_null = findall(isnan.(s.t_infection_active))
    
    # Check expression indices
    @assert all(1 .<= s.expression_index[active_indices] .<= p.n_genes_per_strain)
    @assert all(s.expression_index[active_indices_null] .== 0)
    
    # Check strain IDs
    @assert all(s.strain_id_liver[liver_indices] .> 0)
    @assert all(s.strain_id_liver[liver_indices_null] .== 0)
    @assert all(s.strain_id_active[active_indices] .> 0)
    @assert all(s.strain_id_active[active_indices_null] .== 0)
    
    # Check next strain ID
    @assert s.next_strain_id > max(maximum(s.strain_id_liver), maximum(s.strain_id_active))
    
    # Check genes
    for ci in liver_indices_null
        if !all(s.genes_liver[ci[1], ci[2], :, :] .== 0)
            println("genes_liver[$(ci)] = $(s.genes_liver[ci[1], ci[2], :, :])")
        end
        @assert all(s.genes_liver[ci[1], ci[2], :, :] .== 0)
    end
#     @assert all(s.genes_liver[liver_indices_null, :, :] .== 0)
    for ci in active_indices_null
        @assert all(s.genes_active[ci[1], ci[2], :, :] .== 0)
    end
#     @assert all(s.genes_active[active_indices_null, :, :] .== 0)
    for locus in p.n_loci
        for ci in liver_indices
            @assert all(1 .<= s.genes_liver[ci[1], ci[2], :, locus] .<= s.n_alleles[locus])
        end
        for ci in active_indices
            @assert all(1 .<= s.genes_active[ci[1], ci[2], :, locus] .<= s.n_alleles[locus])
        end
#         @assert all(1 .<= s.genes_liver[liver_indices, :, locus] .<= s.n_alleles[locus])
#         @assert all(1 .<= s.genes_active[active_indices, :, locus] .<= s.n_alleles[locus])
    end
    
    # Check immunity levels
    @assert size(s.immunity)[2] >= maximum(s.n_alleles)
    @assert size(s.immunity)[1] == p.n_hosts
    @assert size(s.immunity)[3] == p.n_loci
    @assert all(s.immunity .< p.immunity_level_max)
    for locus in p.n_loci
        @assert all(s.immunity[:, (s.n_alleles[locus] + 1):size(s.immunity)[2], locus] .== 0)
    end
    
    # Check gene pool
    @assert size(s.gene_pool)[1] .== p.n_genes_initial
    @assert all(s.gene_pool .> 0)
    for locus in 1:p.n_loci
        @assert all(s.gene_pool[:, locus] .<= s.n_alleles[locus])
    end
end

function draw_host_lifetime(p::Params)
    min(rand(Exponential(p.mean_host_lifetime)), p.max_host_lifetime)
end

function gene_similarity(p::Params, gene1, gene2, breakpoint)
    # Directly from C code.
    # TODO: understand it.
    # TODO: factor out these hard-coded constants as parameters.
    p_div = 0.0
    child_div = 0.0
    rho = 0.8
    avg_mutation = 5.0
    for i in 1:p.n_loci
        if gene1[i] != gene2[i]
            p_div += 1.0
            if i < breakpoint
                child_div += 1.0
            end
        end
    end
    rho_power = child_div * avg_mutation * (p_div - child_div) * avg_mutation / (p_div * avg_mutation - 1.0)
    
    rho^rho_power
end

function recombine_genes(gene1, gene2, breakpoint)
    vcat(gene1[1:(breakpoint-1)], gene2[breakpoint:end])
end

function count_infections(s, hosts)
    sum(s.expression_index[hosts, :] .== 0, dims = 2)
end

