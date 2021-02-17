using DataStructures
using Parameters
using Random
using Distributions
using StatsBase

const HostId = UInt16
# const InfectionId = UInt32
const Locus = UInt8
const AlleleId = UInt16 # TODO: 16-bit allele IDs?
const ExpressionIndex = UInt8
const Gene = Array{AlleleId}
const Strain = Array{Union{Nothing, Gene}}
const ImmunityCount = Int8

@with_kw mutable struct Infection
    "Strain: genes in expression order"
    strain::Strain
    
    "Time infection occurred at"
    t::Float32
    
    "Reference to containing Host"
    host_id::HostId # TODO: replace with ID to save memory
    
    "Current expression index of genes"
    expression_index::UInt8
end

const ImmunityCounter = Accumulator{AlleleId, ImmunityCount}

@with_kw mutable struct Host
    "Unique ID across simulation"
    id::HostId
    
    t_birth::Float32
    t_death::Float32
    
    infections::Array{Infection}
    
    "Immunity counts, by locus, by allele"
    immunity_counts::Array{ImmunityCounter}
end

struct ImmunityRef
    host_id::HostId
    locus::Locus
    allele_id::AlleleId
end

@with_kw mutable struct State
    n_alleles::Array{AlleleId}
    gene_pool::Array{Gene}
    
    next_host_id::HostId
    host_id_map::Dict{HostId, Host}
    hosts::IndexedSet{Host}
    infected_hosts::IndexedSet{Host}
    
    infections::IndexedSet{Infection}
    old_infections::Array{Infection}
    
    immunities::IndexedSet{ImmunityRef}
end

struct BitingEventReusableArrays
    inf_indices::Array{Int}
    inf_to_transmit::Array{Int}
    src_strains::Array{Strain}
    strains1::Array{Union{Nothing, Strain}}
    strains2::Array{Strain}
    combined_genes::Strain
    
    function BitingEventReusableArrays()
        new(
            [], [], [], [], [], []
        )
    end
end

function infected_fraction(s::State)
    # NB: In Julia, x / y is always floating-point division, even with integers
    # (x ÷ y or div(x, y) is integer division)
    length(s.infected_hosts) / length(s.hosts)
end

const N_EVENTS = 6
const ALL_EVENTS = collect(1:N_EVENTS)
const (
    BITING,
    IMMIGRATION,
    IMMUNITY_LOSS,
    TRANSITION,
    MUTATION,
    RECOMBINATION
) = ALL_EVENTS

function run_exact(p::Params)
    # println("run_exact()")
    
    rng_seed = if p.rng_seed == nothing
        rng_seed = rand(RandomDevice(), UInt64)
        # println("random seed: $(rng_seed)")
        rng_seed
    else
        p.rng_seed
    end
    Random.seed!(rng_seed)
        
    s = initialize_state(p)
    
    t = 0.0
    t_next_print = 30.0
    
    rates = fill(0.0, N_EVENTS)
    update_rates!(p, t, s, rates, (1:N_EVENTS)...)
    
    exp_dist = Exponential(1.0)
    
    r_biting_event = BitingEventReusableArrays()
    while true
        total_rate = sum(rates)
        dt = rand(exp_dist) / total_rate
        
        t_next = t + dt
        if t_next > p.t_end
            break
        end
        
        t = t_next
        if t > t_next_print
            println("t = $(t)")
            println("frac infected: $(infected_fraction(s))")
            println("# infections: $(length(s.infections))")
            println("# old infections: $(length(s.old_infections))")
            println("# immunities: $(length(s.immunities))")
            t_next_print += 30.0
            
#             Base.GC.gc()
        end
        
        @assert total_rate > 0.0
        
        event_id = sample(1:N_EVENTS, Weights(rates, total_rate))
        if event_id == BITING
            do_biting_event!(p, t, s, r_biting_event)
            update_rates!(
                p, t, s, rates,
                BITING, IMMIGRATION, TRANSITION, MUTATION, RECOMBINATION
            )
        elseif event_id == IMMIGRATION
            do_immigration_event!(p, t, s)
            update_rates!(
                p, t, s, rates,
                BITING, IMMIGRATION, TRANSITION, MUTATION, RECOMBINATION
            )
        elseif event_id == IMMUNITY_LOSS
            do_immunity_loss_event!(p, t, s)
            update_rates!(
                p, t, s, rates,
                BITING, IMMIGRATION
            )
        elseif event_id == TRANSITION
            do_transition_event!(p, t, s)
            update_rates!(
                p, t, s, rates,
                ALL_EVENTS...
            )
        elseif event_id == MUTATION
            do_mutation_event!(p, t, s)
            update_rates!(
                p, t, s, rates,
                BITING, IMMIGRATION
            )
        elseif event_id == RECOMBINATION
            do_recombination_event!(p, t, s)
            update_rates!(
                p, t, s, rates,
                BITING, IMMIGRATION
            )
        end
        
        # println("t = $(t)")
    end
end

function update_rates!(p::Params, t::Float64, s::State, rates, event_ids...)
    for event_id in event_ids
        rates[event_id] = event_rate(event_id, p::Params, t::Float64, s::State)
    end
end

function event_rate(event_id, p::Params, t::Float64, s::State)
    if event_id == BITING
        biting_rate(p, t, s)
    elseif event_id == IMMIGRATION
        immigration_rate(p, t, s)
    elseif event_id == IMMUNITY_LOSS
        immunity_loss_rate(p, t, s)
    elseif event_id == TRANSITION
        transition_rate(p, t, s)
    elseif event_id == MUTATION
        mutation_rate(p, t, s)
    elseif event_id == RECOMBINATION
        recombination_rate(p, t, s)
    end
end

function biting_rate(p::Params, t::Float64, s::State)
    per_capita_rate = p.biting_rate_mean * p.daily_biting_rate_distribution[
        1 + Int(floor(t)) % p.t_year
    ] * infected_fraction(s)
    per_capita_rate * p.n_hosts
end

function do_biting_event!(p::Params, t::Float64, s::State, r::BitingEventReusableArrays)
    # println("do_biting_event!()")
    
    # TODO: is it OK for these to be the same host?
    src_host = rand(s.infected_hosts)
    dst_host = rand(s.hosts)
    
    # Allow hosts to die and get reinitialized if we've passed their death time
    reinitialize_if_past_death!(p, t, s, src_host)
    reinitialize_if_past_death!(p, t, s, dst_host)
    
    if length(dst_host.infections) >= p.max_infection_count
        return
    end
    
    # println("src_host = $(src_host.id), dst_host = $(dst_host.id)")

    # Identify infections past the liver stage
    findall!(r.inf_indices, (inf.t + p.t_liver_stage < t for inf in src_host.infections))
    inf_count = length(r.inf_indices)

    # Only transmit if there are between 1 and moi_transmission_max infections
    if !(1 <= inf_count <= p.moi_transmission_max)
        return
    end
    
    p_infect = p.transmissibility /
        (p.coinfection_reduces_transmission ? 1 : inf_count)
    findall!(r.inf_to_transmit, (rand(inf_count) .< p_infect))
    # println("transmitting $(length(inf_to_transmit)) infections")
    
    collect!(r.src_strains, (inf.strain for inf in src_host.infections[r.inf_to_transmit]))
    
    if length(r.src_strains) == 0
        return
    end
    
    collect!(r.strains1, r.src_strains)
    rand!(r.strains1, r.src_strains)
    
    collect!(r.strains2, r.src_strains)
    rand!(r.strains2, r.src_strains)
    
    for i in 1:min(length(r.src_strains), p.max_infection_count - length(dst_host.infections))
        infection = create_or_reuse_infection!(p, t, s, dst_host)
        
        if r.strains1[i] === r.strains2[i]
            # If they're the same, just shuffle genes
            if length(r.strains1[i]) != length(infection.strain)
                println(r.strains1[i])
                println(infection.strain)
                @assert false
            end
            sample!(r.strains1[i], infection.strain, replace = false)
        else
            # If they're different, sample from both sets of genes without replacement
            vcat!(r.combined_genes, r.strains1[i], r.strains2[i])
            sample!(r.combined_genes, infection.strain, replace = false)
        end
        infect_host!(p, t, s, dst_host, infection)
    end
end

function active_infection_count(p::Params, t::Float64, host::Host)
    sum(infection.t + p.t_liver_stage < t for infection in host.infections)
end


function immigration_rate(p::Params, t::Float64, s::State)
    p.immigration_rate_fraction * biting_rate(p, t, s)
end

function do_immigration_event!(p::Params, t::Float64, s::State)
    @assert p.immigration_on
    host = rand(s.hosts)
    if length(host.infections) < p.max_infection_count
        infection = create_or_reuse_infection!(p, t, s, host)
        rand!(infection.strain, s.gene_pool)
        infect_host!(p, t, s, host, infection)
    end
end

function immunity_loss_rate(p::Params, t::Float64, s::State)
    p.immunity_loss_rate * length(s.immunities)
end

function do_immunity_loss_event!(p::Params, t::Float64, s::State)
    imref = rand(s.immunities)
    host = s.host_id_map[imref.host_id]
    host.immunity_counts[imref.locus][imref.allele_id] -= 1
    
    if host.immunity_counts[imref.locus][imref.allele_id] == 0
        delete!(s.immunities, imref)
    end
end

function transition_rate(p::Params, t::Float64, s::State)
    p.transition_rate * length(s.infections)
end

function do_transition_event!(p::Params, t::Float64, s::State)
    # println("do_transition_event!()")
    
    infection = rand(s.infections)
    host::Host = s.host_id_map[infection.host_id]
    if infection.t + p.t_liver_stage < t
        # println("expression_index = $(infection.expression_index)")
        # Advance expression until non-immune gene is reached
        while true
            # Gain extra immunity to current gene
            gain_immunity!(p, t, s, host, current_gene(infection))
            
            if infection.expression_index == p.n_genes_per_strain
                # Clear infection if we're at the end
                clear_infection!(p, t, s, host, infection)
                break
            else
                # Advance expression if we're not yet at the end
                infection.expression_index += 1
                
                # If we're not immune to this gene, stop advancing expression
                if !is_immune(host, infection)
                    # println("not immune.")
                    break
                end
                # println("immune!")
            end
        end
    end
end

function gain_immunity!(p::Params, t::Float64, s::State, host::Host, gene::Gene)
    # Increment immunity count at each locus
    for locus in 1:lastindex(gene)
        counter = host.immunity_counts[locus]
        allele_id = gene[locus]
        
        if counter[allele_id] < p.max_immunity_count
            inc!(counter, allele_id)
            # Register this host/locus/allele combo for loss by random sampling
            if p.immunity_loss_rate > 0.0 && counter[allele_id] == 1
                push!(s.immunities, ImmunityRef(host.id, locus, allele_id))
            end
        end
    end
end

function current_gene(infection::Infection)
    infection.strain[infection.expression_index]
end

function clear_infection!(
    p::Params, t::Float64, s::State, host::Host, infection::Infection
)
    # println("clear_infection!()")
    delete!(s.infections, infection)
    delete!(host.infections, infection)
    push!(s.old_infections, infection)
    
    if length(host.infections) == 0
        delete!(s.infected_hosts, host)
    end
end

"""
Removes an item in the middle of an array that does not need to be kept ordered
in constant time.

The item is replaced with the item at the end of the array, and then the item
at the end of the array is removed.
"""
function swap_with_end_and_remove!(a, index)
    if index != lastindex(a)
        setindex!(a, a[lastindex(a)], index)
    end
    pop!(a)
    nothing
end

function is_immune(host::Host, infection::Infection)
    gene = current_gene(infection)
    for locus in 1:lastindex(host.immunity_counts)
        allele_id = gene[locus]
        if host.immunity_counts[locus][allele_id] == 0
            return false
        end
    end
    true
end

function mutation_rate(p::Params, t::Float64, s::State)
    p.mutation_rate * p.n_genes_per_strain * p.n_loci * length(s.infections)
end

function do_mutation_event!(p::Params, t::Float64, s::State)
    # println("do_mutation_event!($(t))")
    infection = rand(s.infections)
    
    index = rand(1:p.n_genes_per_strain)
    infection.strain[index] = mutate_gene!(p, t, s, infection.strain[index])
end

function mutate_gene!(p::Params, t::Float64, s::State, gene::Gene)
    gene = Gene(gene)
    locus = rand(1:p.n_loci)
    
    @assert s.n_alleles[locus] < typemax(AlleleId)
    s.n_alleles[locus] += 1
    gene[locus] = s.n_alleles[locus]
    gene
end

function recombination_rate(p::Params, t::Float64, s::State)
    if p.n_genes_per_strain == 1
        0.0
    else
        p.ectopic_recombination_rate * p.n_genes_per_strain *
            (p.n_genes_per_strain - 1) / 2.0 * length(s.infections)
    end
end

function do_recombination_event!(p::Params, t::Float64, s::State)
#     println("do_recombination_event!($(t))")
    
    @assert p.n_genes_per_strain > 1
    
    infection = rand(s.infections)
    
    (index1, index2) = samplepair(p.n_genes_per_strain)
    src_gene_1 = infection.strain[index1]
    src_gene_2 = infection.strain[index2]
    
    if src_gene_1 === src_gene_2
        return
    end
    
    # TODO: need conversions?
    
    breakpoint = rand(1:p.n_loci)
    (new_gene_1, new_gene_2) = if breakpoint == 1
        (src_gene_1, src_gene_2)
    else
        similarity = gene_similarity(p, src_gene_1, src_gene_2, breakpoint)
#         println("similarity = $(similarity)")
        
        new_gene_1 = if rand() < similarity
#             println("1 is functional")
            recombine_genes(src_gene_1, src_gene_2, breakpoint)
        else
            src_gene_1
        end
        
        new_gene_2 = if rand() < similarity
#             println("2 is functional")
            recombine_genes(src_gene_2, src_gene_1, breakpoint)
        else
            src_gene_2
        end
        
        (new_gene_1, new_gene_2)
    end
    
    infection.strain[index1] = new_gene_1
    infection.strain[index2] = new_gene_2
end

function recombine_genes(gene1, gene2, breakpoint)
    vcat(gene1[1:(breakpoint-1)], gene2[breakpoint:end])
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

function initialize_state(p::Params)
    n_alleles = repeat([p.n_alleles_per_locus_initial], p.n_loci)
    gene_pool = initialize_gene_pool(p)
    hosts = initialize_hosts(p, gene_pool)
    s = State(;
        n_alleles = n_alleles,
        gene_pool = gene_pool,
        
        next_host_id = length(hosts) + 1,
        host_id_map = Dict((host.id, host) for host in hosts),
        hosts = hosts,
        infected_hosts = IndexedSet{Host}(),
        
        infections = IndexedSet{Infection}(),
        immunities = IndexedSet{ImmunityRef}(),
        
        old_infections = [],
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
    hosts = IndexedSet{Host}()
    
    for id in 1:p.n_hosts
        lifetime = draw_host_lifetime(p)
        t_birth = -rand() * lifetime
        t_death = t_birth + lifetime
        
        push!(hosts, Host(;
            id = id,
            t_birth = t_birth,
            t_death = t_death,
            infections = Infection[],
            immunity_counts = [ImmunityCounter() for i in 1:p.n_loci]
        ))
    end
    
    hosts
end

function reinitialize_if_past_death!(p::Params, t::Float64, s::State, host::Host)
    while t > host.t_death
#         println("Reinitializing host $(host.id)")
        # Remove all infections
        for infection in host.infections
            delete!(s.infections, infection)
            push!(s.old_infections, infection)
        end
#         host.infections = []
        
        if length(host.infections) > 0
            delete!(s.infected_hosts, host)
        end
        empty!(host.infections)
        
        
        # Remove all immunities
        for locus in 1:p.n_loci
            for (allele_id, count) in host.immunity_counts[locus]
                imref = ImmunityRef(host.id, locus, allele_id)
                if count > 0
                    old_global_immunity_count = length(s.immunities)
                    delete!(s.immunities, imref)
                    @assert length(s.immunities) == old_global_immunity_count - 1
                else
                    @assert !(imref in s.immunities)
                end
            end
#             host.immunity_counts[locus] = ImmunityCounter()
            empty!(host.immunity_counts[locus].map)
        end
        
        # Update ID
        delete!(s.host_id_map, host.id)
        host.id = s.next_host_id
        s.next_host_id += 1
        s.host_id_map[host.id] = host
        
        # Update lifetime
        host.t_birth = host.t_death
        host.t_death = host.t_birth + draw_host_lifetime(p)
    end
end

function draw_host_lifetime(p::Params)
    min(rand(Exponential(p.mean_host_lifetime)), p.max_host_lifetime)
end

function initialize_infections!(p::Params, s::State)
    for id in 1:p.n_initial_infections
        while true
            host = rand(s.hosts)
            
            if length(host.infections) < p.max_infection_count
                infection = create_or_reuse_infection!(p, 0.0, s, host)
                rand!(infection.strain, s.gene_pool)
                infect_host!(p, 0.0, s, host, infection)
                break
            end
        end
    end
end

function infect_host!(p::Params, t::Float64, s::State, host::Host, infection)
    @assert length(host.infections) < p.max_infection_count
    
    # println("Infecting $(host.id)")
    push!(host.infections, infection)
    push!(s.infections, infection)
    
    if length(host.infections) == 1
        push!(s.infected_hosts, host)
    end
end

function create_or_reuse_infection!(p::Params, t::Float64, s::State, host::Host)
    if length(s.old_infections) == 0
        infection = Infection(
            strain = Strain(nothing, p.n_genes_per_strain),
            t = t,
            host_id = host.id,
            expression_index = 1
        )
        
        infection
    else
        infection = pop!(s.old_infections)
        infection.t = t
        infection.host_id = host.id
        infection.expression_index = 1
        
        infection
    end
end
