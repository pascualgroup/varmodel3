"""
Type definitions for model state.

Includes simple integer type definitions for IDs and gene vectors.
If integer limits are exceeded, this code could be modified to use larger
integer types (which could be controlled by a switch in Params).

Also contains `verify()`, which does some cursory checks on state consistency.
"""

const HostId = UInt32
const InfectionId = UInt32
const AlleleId = UInt16
const StrainId = UInt32
const ExpressionIndex = UInt8
const ImmunityLevel = UInt8
const Gene = SVector{P.n_loci, AlleleId}  # Immutable fixed-size vector
const MGene = MVector{P.n_loci, AlleleId} # Mutable fixed-size vector

"""
Struct representing infections.

Infection state currently includes the time at which the infection occurred,
an identifier for the strain (unordered set of genes), the actual genes (as a
matrix of allele IDs), and the currently expressed index.
"""
@with_kw mutable struct Infection
    "Infection identifier, unique across all hosts."
    id::InfectionId
    
    "Time at which infection occurred (entered the liver stage)."
    t_infection::Float64
    
    """
    Strain identifier, unique across the simulation.
    
    When an infection occurs, if the infecting strain does not arise from
    meiotic recombination of two strains, the infecting strain will have the
    same unordered set of genes as one of the strains from the source host,
    but expressed in a new randomly sampled order. In this case, the infection
    is given the same strain ID.
    
    If meiotic recombination occurs, the new strain will be given a new, unique
    strain ID.

    There is a chance that two strains consisting of identical sets of genes
    will arise independently, and thus have different strain IDs.
    Thus, although this representation is more computationally efficient and
    more useful in terms of evolutionary history, if you want to identify
    strains by their genes, you have to actually compare the genes.
    """
    strain_id::StrainId
    
    """
    Genes, specified as an (n_loci, n_genes_per_strain) matrix of `AlleleId`s.
    
    Allele identifiers are unique for a each locus.
    """
    genes::MMatrix{P.n_loci, P.n_genes_per_strain, AlleleId}
    
    """
    Index in `genes` matrix of currently expressed gene.
    
    Set to `0` (and ignored) for liver-stage infections.
    """
    expression_index::ExpressionIndex
end

"""
Struct representing hosts.

Host state currently includes the time at which the host was born; the time at
which the host is scheduled to die; an array of infections currently waiting
in the liver stage; an array of active infections, and immune history.

Infection arrays are dynamically sized but currently limited to
`n_infections_liver_max` / `n_infections_active_max`.
"""
@with_kw mutable struct Host
    "Host identifier, unique across the simulation."
    id::HostId
    
    "Birth time of host."
    t_birth::Float64
    
    "Death time of host."
    t_death::Float64
    
    "Infections in liver stage, not yet activated."
    liver_infections::Array{Infection}
    
    "Actively expressed infections."
    active_infections::Array{Infection}
    
    """
    Immune history.
    
    A dictionary (hash table) mapping genes to immunity level, which are
    stored as 8-bit unsigned integers for memory efficiency and saturate at
    `immunity_level_max`. Expression of a gene results in an incremented
    immunity level; immunity loss results in a decremented immunity level.
    
    When the level reaches 0, the gene is removed from the dictionary.
    """
    immunity::Dict{Gene, ImmunityLevel}
end

"""
Simulation state.

Simulation state includes a global gene pool, used to seed initial infections
and immigration events; an array of hosts; and various accounting and memory
management auxiliaries.
"""
@with_kw mutable struct State
    """
    Gene pool used to generate initial infections and immigration events.
    Entry (i, j) is the allele ID for gene j, locus i.

    Currently fixed-size.
    Dimensions: (n_loci, n_genes_initial)
    """
    gene_pool::Array{AlleleId, 2}
    
    """
    Array of hosts.
    
    Size currently does not change. When a host dies/is reborn, its struct
    is re-used.
    """
    hosts::Array{Host}
    
    "Number of alleles for each locus (epitope)"
    n_alleles::MVector{P.n_loci, AlleleId}

    """
    ID for next host to be born.
    
    Whenever a host is reborn, it is given a new ID.
    """
    next_host_id::HostId
    
    """
    ID for next strain to be created.
    
    Whenever a new strain is created via infection, immigration, mutation,
    or ectopic recombination, it is given the next strain ID, and this
    counter is incremented.
    """
    next_strain_id::StrainId
    
    """
    ID for next infection.

    Whenever a new infection occurs, it is given a new ID.
    """
    next_infection_id::InfectionId
    
    """
    Upper bound on number of immunities per host.
    
    In this code, immunity loss (and everything else) is handled via
    rejection sampling. If every host had seen the same number of genes, we
    could sample a host uniformly at random, and then sample a gene
    uniformly at random.

    But different hosts have seen different numbers of genes. We can still
    sample hosts uniformly at random, but in order to make it so that each
    host/gene combination has the same probability of being chosen, we turn
    the event into a no-op with probability

    `length(host.immunity) / n_immunities_per_host_max`

    where n_immunities_per_host_max is a global upper bound on the number of
    immunities a host can have.
    
    Whenever a host gains immunity, we update the upper bound if the
    host has exceeded it. To prevent unbounded growth of this upper bound,
    we periodically recompute the upper bound by scanning every host.
    The period for these global scans is `P.upper_bound_recomputation_period`.
    """
    n_immunities_per_host_max::Int
    
    """
    Upper bound on number of active infections per host.
    
    See explanation of rejection sampling under `n_immunities_per_host_max`.
    """
    n_active_infections_per_host_max::Int
    
    """
        Array of old infections.
        
        Used to prevent allocation of new infections.
    """
    old_infections::Array{Infection}
end

"""
Verify that various pieces of state are consistent with each other.
"""
function verify(t, s::State)
    println("verify($(t), s)")
    
    for host in s.hosts
        if P.n_infections_liver_max !== missing
            @assert length(host.liver_infections) <= P.n_infections_liver_max
        end
        
        if P.n_infections_active_max !== missing
            @assert length(host.active_infections) <= P.n_infections_active_max
        end
        
        for infection in host.liver_infections
            @assert infection.expression_index == 0
        end
        
        for infection in host.active_infections
            @assert 1 <= infection.expression_index <= P.n_genes_per_strain
        end
    end
end