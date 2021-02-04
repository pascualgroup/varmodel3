using DataStructures
using Parameters
using Random
using Distributions
using StatsBase

const HostId = Int
const HostIndex = Int32
const InfectionId = Int
const InfectionIndex = Int32
const HostInfectionIndex = Int16
const Locus = Int8
const AlleleId = Int
const Gene = Array{AlleleId}
const ImmunityCount = Int8
const ImmunityIndex = Int32

@with_kw mutable struct Infection
    "Time infection occurred"
    t::Float64
    
    "Expression order of genes"
    expression_order::Array{Gene}
    
    "Current expression index of genes"
    expression_index::Int
    
    "Index in global random-samplable array of all infections"
    ref_index::Int
end

struct InfectionRef
    host_index::HostIndex
    infection_index::HostInfectionIndex
end

struct ImmunityRef
    host_index::HostIndex
    locus::Locus
    allele_id::AlleleId
end

@with_kw mutable struct Host
    "Unique ID across simulation"
    id::Int
    
    "Index in global random-samplable array"
    index::HostIndex
    
    t_birth::Float64
    t_death::Float64
    
    infections::Array{Infection}
    
    "Immunity counts, by locus, by allele"
    immunity_counts::Array{Dict{AlleleId, ImmunityCount}}
    
    "Immunity location in global random-samplable array, by locus, by allele"
    immunity_ref_indices::Array{Dict{AlleleId, Int32}}
end

@with_kw mutable struct State
    transmission_count::Int
    bite_count::Int
    infected_bite_count::Int
    
    n_alleles::Array{Int}
    gene_pool::Array{Gene}
    
    hosts::Array{Host}
    
    "Random-samplable array of infected hosts"
    infected_hosts::Array{HostIndex}
    
    "Random-samplable array of all infections in all hosts"
    infections::Array{InfectionRef}
    
    "Random-samplable array of immunity in all hosts"
    immunities::Array{ImmunityRef}
end

function infected_fraction(s::State)
    # NB: In Julia, x / y is always floating-point division, even with integers
    # (x ÷ y or div(x, y) is integer division)
    length(s.infected_hosts) / length(s.hosts)
end

const N_EVENTS = 6
const (
    BITING,
    IMMIGRATION,
    IMMUNITY_LOSS,
    TRANSITION,
    MUTATION,
    RECOMBINATION
) = 1:N_EVENTS

function run_exact(p::Params)
    println("run_exact()")
    
    rng_seed = if p.rng_seed == nothing
        rng_seed = rand(RandomDevice(), UInt64)
        println("random seed: $(rng_seed)")
        rng_seed
    else
        p.rng_seed
    end
    Random.seed!(rng_seed)
        
    s = initialize_state(p)
    
    t = 0.0
    rates = [
        biting_rate(p, t, s),
        immigration_rate(p, t, s),
        immunity_loss_rate(p, t, s),
        transition_rate(p, t, s),
        mutation_rate(p, t, s),
        recombination_rate(p, t, s),
    ]
    
    exp_dist = Exponential(1.0)
    while true
        total_rate = sum(rates)
        dt = rand(exp_dist) * total_rate
        
        t_next = t + dt
        if t_next > p.t_end
            break
        end
        
        t = t_next
        
        @assert total_rate > 0.0
        
        event_id = sample(1:N_EVENTS, Weights(rates, total_rate))
        if event_id == BITING
            do_biting_event!(p, t, s)
        elseif event_id == IMMIGRATION
            do_immigration_event!(p, t, s)
        elseif event_id == IMMUNITY_LOSS
            do_immunity_loss_event!(p, t, s)
        elseif event_id == TRANSITION
            do_transition_event!(p, t, s)
        elseif event_id == MUTATION
            do_mutation_event!(p, t, s)
        elseif event_id == RECOMBINATION
            do_recombination_event!(p, t, s)
        end
        
        println("t = $(t)")
    end
end

function biting_rate(p::Params, t::Float64, s::State)
    per_capita_rate = p.biting_rate_mean * p.daily_biting_rate_distribution[
        1 + Int(floor(t)) % p.t_year
    ] * infected_fraction(s)
    per_capita_rate * p.n_hosts
end

function do_biting_event!(p::Params, t::Float64, s::State)
    println("do_biting_event!()")
    
    # TODO: is it OK for these to be the same host?
    src_host = s.hosts[rand(s.infected_hosts)]
    dst_host = rand(s.hosts)
    
    println("src_host = $(src_host.id), dst_host = $(dst_host.id)")
    
    inf_indices = expression_indices(p, t, src_host)
    active_infections = findall(inf_indices .> 0)
    inf_count = length(inf_indices)
    
    # Only transmit if there are between 1 and moi_transmission_max infections
    if 1 <= inf_count <= p.moi_transmission_max
        p_infect = p.transmissibility / (p.coinfection_reduces_transmission ? 1 : inf_count)
        inf_to_transmit = findall(p_infect .< rand(inf_count))
        println("transmitting $(length(inf_to_transmit)) infections")
        
        src_strains = [inf.expression_order for inf in src_host.infections[inf_to_transmit]]
        
        strains1 = sample(src_strains, length(src_strains))
        strains2 = sample(src_strains, length(src_strains))
        
        for i in 1:length(src_strains)
            strain = if strains1[i] == strains2[i]
                # Shuffle
                sample(strains1[i], p.n_genes_per_strain, replace = false)
            else
                # Random subset of genes in both strains
                sample(vcat(strains1[i], strains2[i]), p.n_genes_per_strain, replace = false)
            end
            infect_host!(t, s, dst_host, strain)
        end
    end
end

function expression_indices(p::Params, t::Float64, host::Host)
    collect(
        inf.t + p.t_liver_stage < t ? inf.expression_index : 0
        for inf in host.infections
    )
end

function active_infection_count(p::Params, t::Float64, host::Host)
    sum(infection.t + p.t_liver_stage < t for infection in host.infections)
end


function immigration_rate(p::Params, t::Float64, s::State)
    p.immigration_rate_fraction * biting_rate(p, t, s)
end

function do_immigration_event!(p::Params, t::Float64, s::State)
    println("do_immigration_event!()")
end

function immunity_loss_rate(p::Params, t::Float64, s::State)
    p.immunity_loss_rate * length(s.immunities)
end

function do_immunity_loss_event!(p::Params, t::Float64, s::State)
    println("do_immunity_loss_event!()")
end

function transition_rate(p::Params, t::Float64, s::State)
    p.transition_rate * length(s.infections)
end

function do_transition_event!(p::Params, t::Float64, s::State)
    println("do_transition_event!()")
end

function mutation_rate(p::Params, t::Float64, s::State)
    p.mutation_rate * p.n_genes_per_strain * p.n_loci * length(s.infections)
end

function do_mutation_event!(p::Params, t::Float64, s::State)
    println("do_mutation_event!()")
end

function recombination_rate(p::Params, t::Float64, s::State)
    p.ectopic_recombination_rate * p.n_genes_per_strain * (p.n_genes_per_strain - 1) / 2.0 * length(s.infections)
end

function do_recombination_event!(p::Params, t::Float64, s::State)
    println("do_recombination_event!()")
end

function initialize_state(p::Params)
    n_alleles = repeat([p.n_alleles_per_locus_initial], p.n_loci)
    gene_pool = initialize_gene_pool(p)
    hosts = initialize_hosts(p, gene_pool)
    s = State(;
        transmission_count = 0,
        bite_count = 0,
        infected_bite_count = 0,
    
        n_alleles = n_alleles,
        gene_pool = gene_pool,
        
        hosts = hosts,
        infected_hosts = [],
        
        infections = [],
        immunities = []
    )
    initialize_infections!(p, s)
    
    s
end

function initialize_gene_pool(p::Params)
    genes_set = Set{Gene}()
    
    # Generate n_genes distinct genes
    for i in 1:p.n_genes_initial
        while true
            gene = rand(1:p.n_alleles_per_locus_initial, p.n_loci)
            if !(gene in genes_set)
                push!(genes_set, gene)
                break
            end
        end
    end
    
    @assert length(genes_set) == p.n_genes_initial
    
    # Return an array
    collect(genes_set)
end

function initialize_hosts(p::Params, gene_pool)
    hosts::Array{Host} = []
    
    for id in 1:p.n_hosts
        lifetime = draw_host_lifetime(p.mean_host_lifetime, p.max_host_lifetime)
        t_birth = -rand() * lifetime
        t_death = t_birth + lifetime
        
        push!(hosts, Host(;
            id = id,
            index = id,
            t_birth = t_birth,
            t_death = t_death,
            infections = [],
            immunity_counts = [Dict() for i in 1:p.n_loci],
            immunity_ref_indices = [Dict() for i in 1:p.n_loci]
        ))
    end
    
    hosts
end

function draw_host_lifetime(mean_lifetime, max_lifetime)
    min(rand(Exponential(mean_lifetime)), max_lifetime)
end

function initialize_infections!(p::Params, s::State)
    for id in 1:p.n_initial_infections
        host = rand(s.hosts)
        strain = rand(s.gene_pool, p.n_genes_per_strain)
        infect_host!(0.0, s, host, strain)
    end
end

function infect_host!(t::Float64, s::State, host::Host, strain::Array{Gene})
    println("Infecting $(host.id)")
    
    ref_index = length(s.infections) + 1
    infection = Infection(;
        t,
        expression_order = strain,
        expression_index = 1,
        ref_index  = ref_index
    )
    push!(host.infections, infection)
    push!(s.infections, InfectionRef(host.index, ref_index))
    
    if length(host.infections) == 1
        push!(s.infected_hosts, host.index)
    end
end
