using DataStructures
using Parameters
using Random
using Distributions
using StatsBase

const HostId = Int64
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
    "Index in global random-samplable array"
    index::Int
    
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
            do_biting_event(p, t, s)
        elseif event_id == IMMIGRATION
            do_immigration_event(p, t, s)
        elseif event_id == IMMUNITY_LOSS
            do_immunity_loss_event(p, t, s)
        elseif event_id == TRANSITION
            do_transition_event(p, t, s)
        elseif event_id == MUTATION
            do_mutation_event(p, t, s)
        elseif event_id == RECOMBINATION
            do_recombination_event(p, t, s)
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

function do_biting_event(p::Params, t::Float64, s::State)
    println("do_biting_event()")
end

function immigration_rate(p::Params, t::Float64, s::State)
    p.immigration_rate_fraction * biting_rate(p, t, s)
end

function do_immigration_event(p::Params, t::Float64, s::State)
    println("do_immigration_event()")
end

function immunity_loss_rate(p::Params, t::Float64, s::State)
    p.immunity_loss_rate * length(s.immunities)
end

function do_immunity_loss_event(p::Params, t::Float64, s::State)
    println("do_immunity_loss_event()")
end

function transition_rate(p::Params, t::Float64, s::State)
    p.transition_rate * length(s.infections)
end

function do_transition_event(p::Params, t::Float64, s::State)
    println("do_transition_event()")
end

function mutation_rate(p::Params, t::Float64, s::State)
    p.mutation_rate * p.n_genes_per_strain * p.n_loci * length(s.infections)
end

function do_mutation_event(p::Params, t::Float64, s::State)
    println("do_mutation_event()")
end

function recombination_rate(p::Params, t::Float64, s::State)
    p.ectopic_recombination_rate * p.n_genes_per_strain * (p.n_genes_per_strain - 1) / 2.0 * length(s.infections)
end

function do_recombination_event(p::Params, t::Float64, s::State)
    println("do_recombination_event()")
end

function initialize_state(p::Params)
    n_alleles = repeat([p.n_alleles_per_locus_initial], p.n_loci)
    gene_pool = initialize_gene_pool(p)
    hosts = initialize_hosts(p, gene_pool)
    (infections, infected_hosts) = initialize_infections(p, hosts, gene_pool)
    
    State(;
        transmission_count = 0,
        bite_count = 0,
        infected_bite_count = 0,
    
        n_alleles = n_alleles,
        gene_pool = gene_pool,
        
        hosts = hosts,
        infected_hosts = infected_hosts,
        
        infections = infections,
        immunities = []
    )
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

function initialize_infections(p::Params, hosts, gene_pool)
    infections::Array{InfectionRef} = []
    infected_hosts::Array{HostIndex} = []
    
    for id in 1:p.n_initial_infections
        host_index = rand(1:p.n_hosts)
        host = hosts[host_index]
        infection = Infection(;
            t = 0.0,
            expression_order = rand(gene_pool, p.n_genes_per_strain),
            expression_index = 0,
            ref_index = id,
        )
        push!(host.infections, infection)
        push!(infections, InfectionRef(host_index, length(host.infections)))
        
        if length(host.infections) == 1
            push!(infected_hosts, host_index)
        end
    end
    
    (infections, infected_hosts)
end
