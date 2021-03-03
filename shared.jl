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
const InfectionCount = Int8

@with_kw mutable struct State
    # Number of alleles at each locus
    n_alleles::Vector{AlleleId}
    
    t_birth::Vector{Float32}
    t_death::Vector{Float32}
    
    # Current expression index of each active infection
    # host X infection
    expression_index::Array{ExpressionIndex, 2}
    
#     n_infections_liver::Array{InfectionCount}
#     n_infections_active::Array{InfectionCount}
    
    # Time of each infection
    # host X infection
    t_infection_liver::Array{Float32, 2}
    t_infection_active::Array{Float32, 2}
    
    # Strain IDs for each infection, used during infection to
    # compare whether two sets of genes are copies of each other
    # host X infection
    infection_strain_id_liver::Array{StrainId, 2}
    infection_strain_id_active::Array{StrainId, 2}
    
    next_strain_id::StrainId
    
    # Genes that make up each infection, stored directly as
    # sequences of sequences of allele IDs
    # host X infection X expression_index X locus
    infection_genes_liver::Array{AlleleId, 4}
    infection_genes_active::Array{AlleleId, 4}
    
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
    
#     n_infections_liver = fill(0, p.n_hosts)
#     n_infections_active = fill(0, p.n_hosts)
    
    t_infection_liver = fill(NaN32, p.n_hosts, p.infection_count_liver_max)
    t_infection_active = fill(NaN32, p.n_hosts, p.infection_count_active_max)
    infection_strain_id_liver = fill(0, p.n_hosts, p.infection_count_liver_max)
    infection_strain_id_active = fill(0, p.n_hosts, p.infection_count_active_max)
    infection_genes_liver = fill(AlleleId(0), p.n_hosts, p.infection_count_liver_max, p.n_genes_per_strain, p.n_loci)
    infection_genes_active = fill(AlleleId(0), p.n_hosts, p.infection_count_active_max, p.n_genes_per_strain, p.n_loci)
    expression_index = fill(ExpressionIndex(0), p.n_hosts, p.infection_count_active_max)
    
    gene_pool = reshape(rand(
        AlleleId(1):AlleleId(p.n_alleles_per_locus_initial),
        p.n_genes_initial * p.n_loci
    ), (p.n_genes_initial, p.n_loci))
    
    infection_hosts = sample(1:p.n_hosts, p.n_initial_infections, replace = false)
#     n_infections[infection_hosts] .= 1
#     n_infections_liver[infection_hosts] .= 1
    t_infection_liver[infection_hosts] .= 0.0
    
    infection_strain_id_liver[infection_hosts, 1] = 1:p.n_initial_infections
    next_strain_id = p.n_initial_infections + 1
    
    for i in 1:p.n_genes_per_strain
        infection_genes_liver[infection_hosts, 1, i, :] = gene_pool[rand(1:p.n_genes_initial, p.n_initial_infections), :]
    end
    
    State(
        t_birth = t_birth,
        t_death = t_death,
        
#         n_infections = n_infections,
        t_infection_liver = t_infection_liver,
        t_infection_active = t_infection_active,
        
        infection_strain_id_liver = infection_strain_id_liver,
        infection_strain_id_active = infection_strain_id_active,
        
        next_strain_id = next_strain_id,
        
        infection_genes_liver = infection_genes_liver,
        infection_genes_active = infection_genes_active,
        
        expression_index = expression_index,
        
        n_alleles = fill(AlleleId(p.n_alleles_per_locus_initial), p.n_loci),
        immunity = fill(ImmunityCount(0), p.n_hosts, p.n_alleles_per_locus_initial, p.n_loci),
        
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
    @assert all(s.infection_strain_id_liver[liver_indices] .> 0)
    @assert all(s.infection_strain_id_liver[liver_indices_null] .== 0)
    @assert all(s.infection_strain_id_active[active_indices] .> 0)
    @assert all(s.infection_strain_id_active[active_indices_null] .== 0)
    
    # Check next strain ID
    @assert s.next_strain_id > max(maximum(s.infection_strain_id_liver), maximum(s.infection_strain_id_active))
    
    # Check genes
    @assert all(s.infection_genes_liver[liver_indices_null, :, :] .== 0)
    @assert all(s.infection_genes_active[active_indices_null, :, :] .== 0)
    for locus in p.n_loci
        @assert all(1 .<= s.infection_genes_liver[liver_indices, locus, :] .<= s.n_alleles[locus])
        @assert all(1 .<= s.infection_genes_active[active_indices, locus, :] .<= s.n_alleles[locus])
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

function count_infections(s, hosts)
    sum(s.expression_index[hosts, :] .== 0, dims = 2)
end

function mask_with_row_limits(mask, limits)
    count = fill(0, length(limits))
    new_mask = falses(size(mask))
    for i in size(mask)[2]
        new_mask[:,i] = mask[:,i] .& (count .< limits)
        count[:] .+= new_mask[:,i]
    end
    new_mask
end
