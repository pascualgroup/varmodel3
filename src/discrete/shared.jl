using DataStructures
using Parameters
using Random
using Distributions
using StatsBase
using StaticArrays
using Dates
using Test

const HostId = UInt32
const InfectionId = UInt32
const HostGeneIndex = UInt16
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
        Entry (i, j) is the allele ID for gene j, locus i.

        Currently fixed-size.
        Dimensions: (n_loci, n_genes_initial)
    """
    gene_pool::Array{AlleleId, 2}

    """
        ID for next host to be born.
        
        Whenever a host is reborn, it is given a new ID.
    """
    next_host_id::HostId
    
    """
        ID for next strain to be created.

        Whenever a new strain is created via infection, immigration, mutation,
        or ectopic recombination, it is given the next strain ID and this
        counter is incremented.
    """
    next_strain_id::StrainId
    
    """
        ID for next infection.

        Whenever a new infection occurs, it's given a new ID.
    """
    next_infection_id::InfectionId
    
    """
        ID for each host.
        
        Dimensions: (n_hosts,)
    """
    host_id::Vector{HostId}
    
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
        Infection ID for infections in liver stage.

        Dimensions: (n_infections_liver_max, n_hosts)
    """
    infection_id_liver::Array{InfectionId, 2}

    """
        Infection ID for active infections.

        Dimensions: (n_infections_active_max, n_hosts)
    """
    infection_id_active::Array{InfectionId, 2}

    """
        Time of infection for infections in the liver stage.

        Fixed-size array for the sake of predictable memory usage.
        Unused slots are indicated by NaN.

        Dimensions: (n_infections_liver_max, n_hosts)
    """
    t_infection_liver::Array{Float32, 2}

    """
        Time of infection for active infections.

        Fixed-size array for the sake of predictable memory usage.
        Unused slots are indicated by NaN.

        Dimensions: (n_infections_active_max, n_hosts)
    """
    t_infection_active::Array{Float32, 2}

    """
        Strain IDs for infections in the liver stage.

        Whenever a strain is created or modified, it is given a new strain ID.
        Unused slots are indicated by 0.

        Dimensions: (n_infections_liver_max, n_hosts)
    """
    strain_id_liver::Array{StrainId, 2}

    """
        Strain IDs for active infections, transferred directly from the liver stage.

        Unused slots are indicated by 0.

        Dimensions: (n_infections_active_max, n_hosts)
    """
    strain_id_active::Array{StrainId, 2}
    
    """
        Genes for infections in the liver stage, as raw allele IDs.

        Dimensions: (n_loci, n_genes_per_strain, n_infections_liver_max, n_hosts)
    """
    genes_liver::Array{HostGeneIndex, 4}

    """
        Genes for active infections, as raw allele IDs

        Dimensions: (n_loci, n_genes_per_strain, n_infections_active_max, n_hosts)
    """
    genes_active::Array{HostGeneIndex, 4}

    """
        Index of currently expressed gene for each active infection.

        Dimensions: (n_infections_active_max, n_hosts)
    """
    expression_index::Array{ExpressionIndex, 2}
    
    """
        Immunity for each host, as dict mapping alleles to immunity level.
    """
    immunity::Array{Dict}
end

function State(p::Params)
    lifetime = Vector([draw_host_lifetime(p) for i in 1:p.n_hosts])
    t_birth = -rand(Float32, p.n_hosts) .* lifetime
    t_death = t_birth + lifetime
    
    host_id = 1:p.n_hosts
    next_host_id = p.n_hosts + 1
    
    infection_id_liver = fill(InfectionId(0), p.n_infections_liver_max, p.n_hosts)
    infection_id_active = fill(InfectionId(0), p.n_infections_active_max, p.n_hosts)
    t_infection_liver = fill(NaN32, p.n_infections_liver_max, p.n_hosts)
    t_infection_active = fill(NaN32, p.n_infections_active_max, p.n_hosts)
    strain_id_liver = fill(StrainId(0), p.n_infections_liver_max, p.n_hosts)
    strain_id_active = fill(StrainId(0), p.n_infections_active_max, p.n_hosts)
    genes_liver = fill(AlleleId(0), p.n_loci, p.n_genes_per_strain, p.n_infections_liver_max, p.n_hosts)
    genes_active = fill(AlleleId(0), p.n_loci, p.n_genes_per_strain, p.n_infections_liver_max, p.n_hosts)
    expression_index = fill(ExpressionIndex(0), p.n_infections_active_max, p.n_hosts)
    immunity = [Dict{SVector{p.n_loci, AlleleId}, ImmunityLevel}() for i in 1:p.n_hosts]
    
    gene_pool = reshape(rand(
        AlleleId(1):AlleleId(p.n_alleles_per_locus_initial),
        p.n_genes_initial * p.n_loci
    ), (p.n_loci, p.n_genes_initial))
    
    infection_hosts = sample(1:p.n_hosts, p.n_initial_infections, replace = false)
    t_infection_liver[1, infection_hosts] .= 0.0
    infection_id_liver[1, infection_hosts] = 1:p.n_initial_infections
    next_infection_id = p.n_initial_infections + 1

    strain_id_liver[1, infection_hosts] = 1:p.n_initial_infections
    next_strain_id = p.n_initial_infections + 1
    
    genes_liver[:, :, 1, infection_hosts] = reshape(
        gene_pool[:, rand(1:p.n_genes_initial, p.n_genes_per_strain * p.n_initial_infections)],
        (p.n_loci, p.n_genes_per_strain, p.n_initial_infections)
    )

    State(
        next_strain_id = next_strain_id,
        next_host_id = next_host_id,
        next_infection_id = next_infection_id,
        
        host_id = host_id,
        t_birth = t_birth,
        t_death = t_death,
        
        infection_id_liver = infection_id_liver,
        infection_id_active = infection_id_active,
        
        t_infection_liver = t_infection_liver,
        t_infection_active = t_infection_active,

        strain_id_liver = strain_id_liver,
        strain_id_active = strain_id_active,

        genes_liver = genes_liver,
        genes_active = genes_active,

        expression_index = expression_index,
        
        n_alleles = fill(AlleleId(p.n_alleles_per_locus_initial), p.n_loci),
        
        gene_pool = gene_pool,
        
        immunity = immunity
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
    
    # Check infection IDs
    @assert all(s.infection_id_liver[liver_indices] .!= 0)
    @assert all(s.infection_id_liver[liver_indices_null] .== 0)
    @assert all(s.infection_id_active[active_indices] .!= 0)
    @assert all(s.infection_id_active[active_indices_null] .== 0)
    
#     println("liver_indices = $(liver_indices)")

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
    
    # Check references to host genes
    for i in 1:p.n_loci
        @assert all(1 .<= s.genes_liver[i, :, liver_indices] .<= s.n_alleles[i])
        @assert all(1 .<= s.genes_active[i, :, active_indices] .<= s.n_alleles[i])
    end
    @assert all(s.genes_liver[:, :, liver_indices_null] .== 0)
    @assert all(s.genes_active[:, :, active_indices_null] .== 0)
    
    for i in 1:p.n_hosts
        for (gene, level) in s.immunity[i]
            @assert level > 0
        end
    end

    # Check gene pool
    @assert size(s.gene_pool)[2] .== p.n_genes_initial
    @assert all(s.gene_pool .> 0)
    for i in 1:p.n_loci
        @assert all(s.gene_pool[i, :] .<= s.n_alleles[i])
    end
end

function next_strain_id!(s)
    strain_id = s.next_strain_id
    s.next_strain_id += 1
    strain_id
end

function next_strain_ids!(s, n)
    strain_ids = s.next_strain_id:StrainId(s.next_strain_id + n - 1)
    s.next_strain_id += n
    strain_ids
end

function next_host_id!(s)
    next_host_id = s.next_host_id
    s.next_host_id += 1
    host_id
end

function next_host_ids!(s, n)
    host_ids = s.next_host_id:HostId(s.next_host_id + n - 1)
    s.next_host_id += n
    host_ids
end

function next_infection_id!(s)
    infection_id = s.next_infection_id
    s.next_infection_id += 1
    infection_id
end

function next_infection_ids!(s, n)
    infection_ids = s.next_infection_id:InfectionId(s.next_infection_id + n - 1)
    s.next_infection_id += n
    infection_ids
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

function n_host_genes_max(s::State)
    size(s.host_genes)[2]
end

function cartesian_indices_infection_host(s::State)
    CartesianIndices(s.expression_index)
end