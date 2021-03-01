using DataStructures
using Parameters
using Random
using Distributions
using StatsBase
using StaticArrays

const StrainId = UInt32
const AlleleId = UInt16
const Locus = UInt8
const ExpressionIndex = UInt8
const ImmunityCount = Int8

@with_kw mutable struct State
    # Number of alleles at each locus
    n_alleles::Vector{AlleleId}
    
    t_birth::Vector{Float32}
    t_death::Vector{Float32}
    
    # Current expression index of each infection
    # host X infection
    expression_index::Array{ExpressionIndex, 2}
    
    # Time of each infection
    # host X infection
    t_infection::Array{Float32, 2}
    
    # Strain IDs for each infection, used during infection to
    # compare whether two sets of genes are copies of each other
    # host X infection
    infection_strain_id::Array{StrainId, 2}
    
    next_strain_id::StrainId
    
    # Genes that make up each infection, stored directly as
    # sequences of sequences of allele IDs
    # host X infection X expression_index X locus
    infection_genes::Array{AlleleId, 4}
    
    # Immunity level for every allele
    # host X allele_id X locus
    immunity::Array{ImmunityCount, 3}
    
    # Gene pool for immigration events
    # gene X locus
    gene_pool::Array{AlleleId, 2}
end

function State(p::Params)
    lifetime = Vector([draw_host_lifetime(p) for i in 1:p.n_hosts])
    t_birth = -rand(Float32, p.n_hosts) .* lifetime
    t_death = t_birth + lifetime
    
#     n_infections = fill(UInt8(0), p.n_hosts)
    t_infection = fill(NaN32, p.n_hosts, p.max_infection_count)
    infection_strain_id = fill(0, p.n_hosts, p.max_infection_count)
    infection_genes = fill(AlleleId(0), p.n_hosts, p.max_infection_count, p.n_genes_per_strain, p.n_loci)
    expression_index = fill(ExpressionIndex(0), p.n_hosts, p.max_infection_count)
    
    gene_pool = reshape(rand(
        AlleleId(1):AlleleId(p.n_alleles_per_locus_initial),
        p.n_genes_initial * p.n_loci
    ), (p.n_genes_initial, p.n_loci))
    
    infection_hosts = sample(1:p.n_hosts, p.n_initial_infections, replace = false)
#     n_infections[infection_hosts] .= 1
    t_infection[infection_hosts] .= 0.0
    
    infection_strain_id[infection_hosts, 1] = 1:p.n_initial_infections
    next_strain_id = p.n_initial_infections + 1
    
    for i in 1:p.n_genes_per_strain
        infection_genes[infection_hosts, 1, i, :] = gene_pool[rand(1:p.n_genes_initial, p.n_initial_infections), :]
    end
    expression_index[infection_hosts, 1] .= 1
    
    State(
        t_birth = t_birth,
        t_death = t_death,
        
#         n_infections = n_infections,
        t_infection = t_infection,
        
        infection_strain_id = infection_strain_id,
        next_strain_id = next_strain_id,
        
        infection_genes = infection_genes,
        
        expression_index = expression_index,
        
        n_alleles = fill(AlleleId(p.n_alleles_per_locus_initial), p.n_loci),
        immunity = fill(ImmunityCount(0), p.n_hosts, p.n_alleles_per_locus_initial, p.n_loci),
        
        gene_pool = gene_pool
    )
end

function verify(p::Params, s::State)
    @assert all(.!isnan.(s.t_birth))
    @assert all(.!isnan.(s.t_death))
    
    # Check expression indices; use to identify active and inactive infections
    @assert all(s.expression_index .>= 0 .| s.expression_index .<= p.n_genes_per_strain)
    inactive_inf_indices = findall(s.expression_index .== 0)
    active_inf_indices = findall(s.expression_index .!= 0)
    
    # Check that infection time NaNs match expression index 0's
    @assert all(isnan.(s.t_infection[inactive_inf_indices]))
    @assert all(.!isnan.(s.t_infection[active_inf_indices]))
    
    # Check that strain ID 0's match expression index 0's
    @assert all(s.infection_strain_id[inactive_inf_indices] .== 0)
    @assert all(s.infection_strain_id[active_inf_indices] .> 0)
    
    # Check next strain ID
    @assert s.next_strain_id > maximum(s.infection_strain_id)
    
    # Check active and inactive infection genes
    @assert all(s.infection_genes[inactive_inf_indices, :, :] .== 0)
    @assert all(s.infection_genes[active_inf_indices, :, :] .> 0)
    for locus in p.n_loci
        @assert all(s.infection_genes[active_inf_indices, locus, :] .<= s.n_alleles[locus])
    end
    
    # Check immunity levels
    @assert size(s.immunity)[2] >= maximum(s.n_alleles)
    @assert size(s.immunity)[1] == p.n_hosts
    @assert size(s.immunity)[3] == p.n_loci
    @assert all(s.immunity .< p.max_immunity_count)
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
