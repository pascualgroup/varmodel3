const HostId = UInt32
const InfectionId = UInt32
const AlleleId = UInt16
const StrainId = UInt32
const ExpressionIndex = UInt8
const ImmunityLevel = UInt8

@with_kw mutable struct Infection
    id::InfectionId
    t_infection::Float64
    strain_id::StrainId
    genes::Matrix{AlleleId}
    expression_index::ExpressionIndex
end

@with_kw mutable struct Host
    id::HostId
    t_birth::Float64
    t_death::Float64
    liver_infections::Array{Infection}
    active_infections::Array{Infection}
    immunity::Dict{SVector{P.n_loci, AlleleId}, ImmunityLevel}
end

@with_kw mutable struct State
    "Number of alleles for each locus (epitope)"
    n_alleles::MVector{P.n_loci, AlleleId}

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
        Array of hosts.
        
        Size currently does not change. When a host dies/is reborn, its struct
        is re-used.
    """
    hosts::Array{Host}
    
    """
        Upper bound on number of immunities per host.
        
        Updated periodically to prevent excessive growth.
    """
    n_immunities_per_host_max::Int
    
    """
        Array of old infections.
        
        Used to prevent allocation of new infections.
    """
    old_infections::Array{Infection}
end


function verify(t, s)
    println("verify($(t), s)")
    
    for host in s.hosts
        if length(host.liver_infections) > P.n_infections_liver_max
            println("host = $(host.id), n_liver = $(length(host.liver_infections))")
        end
        
        @assert length(host.liver_infections) <= P.n_infections_liver_max
        @assert length(host.active_infections) <= P.n_infections_active_max
        
        for infection in host.liver_infections
            @assert infection.expression_index == 0
        end
        
        for infection in host.active_infections
            @assert 1 <= infection.expression_index <= P.n_genes_per_strain
        end
    end
end
