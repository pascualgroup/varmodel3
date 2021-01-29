using Parameters
using Random
using Distributions

const AlleleId = Int
const Gene = Array{AlleleId}
const ImmunityCounter = Dict{AlleleId, Int}

const Strain = Array{Gene}

@with_kw mutable struct Infection
    t::Float64
    expression_order::Strain
    expression_index::Int
end

@with_kw mutable struct Host
    t_birth::Float64
    t_death::Float64
    
    infections::Array{Infection}
    immunity::Array{ImmunityCounter}
end

@with_kw mutable struct State
    transmission_count::Int
    bite_count::Int
    infected_bite_count::Int
    
    n_alleles::Array{Int}
    gene_pool::Array{Gene}
    hosts::Array{Host}
    infections::Array{Infection}
end

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
        
    state = initialize_state(p)
end

function initialize_state(p::Params)
    n_alleles = repeat([p.n_alleles_per_locus_initial], p.n_loci)
    gene_pool = initialize_gene_pool(p)
    hosts = initialize_hosts(p, gene_pool)
    infections = initialize_infections(p, hosts, gene_pool)
    
    State(;
        transmission_count = 0,
        bite_count = 0,
        infected_bite_count = 0,
    
        n_alleles = n_alleles,
        gene_pool = gene_pool,
        hosts = hosts,
        infections = infections,
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
            t_birth = t_birth,
            t_death = t_death,
            infections = [],
            immunity = [ImmunityCounter() for i in 1:p.n_loci]
        ))
    end
    
    hosts
end

function draw_host_lifetime(mean_lifetime, max_lifetime)
    min(rand(Exponential(mean_lifetime)), max_lifetime)
end

function initialize_infections(p::Params, hosts, gene_pool)
    infections::Array{Infection} = []
    
    for id in 1:p.n_initial_infections
        host = rand(hosts)
        infection = Infection(;
            t = 0.0,
            expression_order = rand(gene_pool, p.n_genes_per_strain),
            expression_index = 0,
        )
        push!(host.infections, infection)
        push!(infections, infection)
    end
    
    infections
end
