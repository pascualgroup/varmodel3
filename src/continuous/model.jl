include("state.jl")
include("output.jl")

const N_EVENTS = 6
const EVENTS = collect(1:N_EVENTS)
const (BITING, IMMIGRATION, SWITCHING, MUTATION, ECTOPIC_RECOMBINATION, IMMUNITY_LOSS) = EVENTS

function run()
    db = initialize_database()
    
    # Seed the random number generator using the provided seed,
    # or, if absent, by generating one from the OS's source of entropy.
    rng_seed = if P.rng_seed === missing
        rand(RandomDevice(), 1:typemax(Int64))
    else
        P.rng_seed
    end
    Random.seed!(rng_seed)
    execute(db.meta, ("rng_seed", rng_seed))
    
    # Start recording elapsed time
    start_datetime = now()
    
    # Used to decide when to do output and verification
    t_next_integer = 1
    
    # Initialize state
    t = 0.0
    s = initialize_state()
    verify(t, s)
    
    # Run initial output
    stats = SummaryStats()
    write_output!(db, 0, s, stats)
    
    # Initialize event rates
    rates = [get_rate(t, s, event) for event in EVENTS]
    total_rate = sum(rates)
    
    # Loop events until end of simulation
    while total_rate > 0.0 && t < P.t_end
        # Draw next time with rate equal to the sum of all event rates
        dt = rand(Exponential(1.0 / total_rate))
        @assert dt > 0.0 && !isinf(dt)
        
        # At each integer time, write output/state verification (if necessary),
        # and update the biting/immigration rate.
        # Loop required in case the simulation jumps past two integer times.
        while t_next_integer < t + dt
            write_output!(db, t_next_integer, s, stats)
            if P.verification_period != nothing && t_next_integer % P.verification_period == 0
                verify(t_next_integer, s)
            end
            
            if t_next_integer % P.upper_bound_recomputation_period == 0
                recompute_rejection_upper_bounds!(s)
            end
            
            t_next_integer += 1
        end
        
        # Draw the event, update time, and execute event
        event = direct_sample_linear_scan(rates, total_rate)
        t += dt
        event_happened = do_event!(t, s, stats, event)
        
        # Many events are no-ops due to rejection sampling.
        # If an event happened, update all rates.
        # (This is wasteful but not a bottleneck.)
        if event_happened
            total_rate = update_rates!(rates, t, s)
            stats.n_events += 1
        end
    end
    
    elapsed_time = Dates.value(now() - start_datetime) / 1000.0
    println("elapsed time (s): $(elapsed_time)")
    execute(db.meta, ("elapsed_time", elapsed_time))
    
    went_extinct = total_rate == 0.0
    println("went extinct? $(went_extinct)")
    execute(db.meta, ("went_extinct", Int64(went_extinct)))
end

function recompute_rejection_upper_bounds!(s)
    s.n_immunities_per_host_max = maximum(length(host.immunity) for host in s.hosts)
    s.n_active_infections_per_host_max = maximum(length(host.active_infections) for host in s.hosts)
end


### INITIALIZATION ###

function initialize_state()
    println("initialize_state()")
    
    # Initialize gene pool as an (n_loci, n_genes_initial) matrix filled with
    # allele IDs drawn uniformly randomly in 1:n_alleles_per_locus_initial.
    gene_pool = reshape(
        rand(1:P.n_alleles_per_locus_initial, P.n_loci * P.n_genes_initial),
        (P.n_loci, P.n_genes_initial)
    )
    
    # Initialize n_hosts hosts, all born at t = 0, with lifetime drawn from a
    # distribution, and no initial infections or immunity.
    hosts = [
        Host(
            id = id,
            t_birth = 0.0, t_death = draw_host_lifetime(),
            liver_infections = [], active_infections = [], immunity = Dict()
        )
        for id in 1:P.n_hosts
    ]
    
    # Infect n_initial_infections hosts at t = 0. Genes in the infection
    # are sampled uniformly randomly from the gene pool.
    for (infection_id, host_index) in enumerate(sample(1:P.n_hosts, P.n_initial_infections, replace = false))
        infection = create_empty_infection()
        infection.id = infection_id
        infection.t_infection = 0.0
        infection.strain_id = infection_id
        infection.genes[:,:] = reshape(
            gene_pool[:, rand(1:P.n_genes_initial, P.n_genes_per_strain)],
            (P.n_loci, P.n_genes_per_strain)
        )
        push!(hosts[host_index].liver_infections, infection)
    end
    
    State(
        n_alleles = fill(P.n_alleles_per_locus_initial, P.n_loci),
        gene_pool = gene_pool,
        next_host_id = P.n_hosts + 1,
        next_strain_id = P.n_initial_infections + 1,
        next_infection_id = P.n_initial_infections + 1,
        hosts = hosts,
        old_infections = [],
        n_immunities_per_host_max = 0,
        n_active_infections_per_host_max = 0
    )
end

"""
    Draws host lifetime from a distribution.
    
    The distribution is an exponential distribution with mean
    `mean_host_lifetime`, truncated at `max_host_lifetime`.
"""
function draw_host_lifetime()
    dist = Exponential(P.mean_host_lifetime)
    while true
        lifetime = rand(dist)
        if lifetime < P.max_host_lifetime
            return lifetime
        end
    end
end

function create_empty_infection()
    Infection(
        id = 0,
        t_infection = NaN,
        strain_id = StrainId(0),
        genes = fill(AlleleId(0), (P.n_loci, P.n_genes_per_strain)),
        expression_index = ExpressionIndex(0)
    )
end

function recycle_or_create_infection(s::State)
    if !isempty(s.old_infections)
        pop!(s.old_infections)
    else
        create_empty_infection()
    end
end


### EVENT DEMUX ###

function update_rates!(rates, t, s)
    for event in 1:N_EVENTS
        rates[event] = get_rate(t, s, event)
    end
    sum(rates)
end

function get_rate(t, s, event)
    if event == BITING
        get_rate_biting(t, s)
    elseif event == IMMIGRATION
        get_rate_immigration(t, s)
    elseif event == SWITCHING
        get_rate_switching(t, s)
    elseif event == MUTATION
        get_rate_mutation(t, s)
    elseif event == ECTOPIC_RECOMBINATION
        get_rate_ectopic_recombination(t, s)
    elseif event == IMMUNITY_LOSS
        get_rate_immunity_loss(t, s)
    end
end

function do_event!(t, s, stats, event)
    if event == BITING
        do_biting!(t, s, stats)
    elseif event == IMMIGRATION
        do_immigration!(t, s, stats)
    elseif event == SWITCHING
        do_switching!(t, s, stats)
    elseif event == MUTATION
        do_mutation!(t, s, stats)
    elseif event == ECTOPIC_RECOMBINATION
        do_ectopic_recombination!(t, s, stats)
    elseif event == IMMUNITY_LOSS
        do_immunity_loss!(t, s, stats)
    end
end


### BITING EVENT ###

function get_rate_biting(t, s)
    biting_rate = P.biting_rate[1 + Int(floor(t)) % P.t_year]
    biting_rate * P.n_hosts
end

function do_biting!(t, s, stats)
#     println("do_biting!($(t), s)")
    
    stats.n_bites += 1
    
    # Uniformly randomly sample infecting host (source) and host being infected
    # (destination).
    src_host = rand(s.hosts)
    dst_host = rand(s.hosts)
    
    # Advance host (rebirth or infection activation)
    advance_host!(t, s, src_host)
    advance_host!(t, s, dst_host)
    
    # The source host must be infected in order to transmit.
    src_active_count = length(src_host.active_infections)
    if src_active_count == 0
        return false
    end
    stats.n_infected_bites += 1
    
    # The destination host must have space available in the liver stage.
    dst_available_count = if P.n_infections_liver_max === missing
        src_active_count
    else
        P.n_infections_liver_max - length(dst_host.liver_infections)
    end
    if dst_available_count == 0
        return false
    end
    stats.n_infected_bites_with_space += 1
    
    # Compute probability of each transmission
    p_transmit = if P.coinfection_reduces_transmission
        P.transmissibility / src_active_count
    else
        P.transmissibility
    end
    
    # The number of transmissions is bounded by the number of source infections
    # and the number of available slots in the destination.
    n_transmissions_max = min(src_active_count, dst_available_count)
    transmitted = false
    for i in 1:n_transmissions_max
        if rand() < p_transmit
            stats.n_transmissions += 1
            transmitted = true
            
            # Randomly sample two source infections to recombine
            src_inf_1 = rand(src_host.active_infections)
            src_inf_2 = rand(src_host.active_infections)
        
            # Get a new infection struct, or recycle an old infection
            # to prevent excess memory allocation.
            dst_inf = recycle_or_create_infection(s)
            dst_inf.id = next_infection_id!(s)
            dst_inf.t_infection = t
            dst_inf.expression_index = 0
            
            # Construct strain for new infection
            if src_inf_1.strain_id == src_inf_2.strain_id
                # If both infections have the same strain, then the new infection
                # is given infection 1's genes with expression order shuffled.
                dst_inf.strain_id = src_inf_1.strain_id
                shuffle_columns_to!(dst_inf.genes, src_inf_1.genes)
            else
                # Otherwise, the new infection is given a new strain constructed by
                # taking a random sample of the genes in the two source infections.
                dst_inf.strain_id = next_strain_id!(s)
                sample_columns_from_two_matrices_to!(dst_inf.genes, src_inf_1.genes, src_inf_2.genes)
            end
            
            # Add this infection to the destination host
            push!(dst_host.liver_infections, dst_inf)
        end
    end
    
    if transmitted
        stats.n_transmitting_bites += 1
        true
    else
        false
    end
end

function advance_host!(t, s, host)
    if t > host.t_death
        # If the host is past its death time
        do_rebirth!(t, s, host)
    else
        i = 1
        while i <= length(host.liver_infections)
            infection = host.liver_infections[i]
            if infection.t_infection + P.t_liver_stage < t
#                 println("t = $(t): activating host $(host.id), inf $(infection.id)")
                
                # If the infection is past the liver stage, remove it from the liver.
                delete_and_swap_with_end!(host.liver_infections, i)
                # If there's room, move it into the active infections array.
                # Otherwise, just put it into the recycle bin.
                if P.n_infections_active_max === missing || length(host.active_infections) < P.n_infections_active_max
                    infection.expression_index = 1
                    push!(host.active_infections, infection)
                    
                    if length(host.active_infections) > s.n_active_infections_per_host_max
                        s.n_active_infections_per_host_max = length(host.active_infections)
                    end
                else
                    push!(s.old_infections, infection)
                end
            else
                i += 1
            end
        end
    end
    nothing
end

function do_rebirth!(t, s, host)
    host.id = next_host_id!(s)
    host.t_birth = t
    host.t_death = t + draw_host_lifetime()
    empty!(host.liver_infections)
    empty!(host.active_infections)
    empty!(host.immunity)
end


### IMMIGRATION EVENT ###

function get_rate_immigration(t, s)
    P.immigration_rate_fraction * get_rate_biting(t, s)
end

function do_immigration!(t, s, stats)
#     println("do_immigration!($(t))")
    
    # Sample a random host and advance it (rebirth or infection activation)
    host = rand(s.hosts)
    advance_host!(t, s, host)
    
    # If host doesn't have an available infection slot, reject this sample.
    if P.n_infections_liver_max !== missing
        if length(host.liver_infections) == P.n_infections_liver_max
            return false
        end
    end
    
    # Construct infection by sampling from gene pool
    infection = recycle_or_create_infection(s)
    infection.id = next_infection_id!(s)
    infection.t_infection = t
    infection.strain_id = next_strain_id!(s)
    infection.expression_index = 0
    for i in 1:P.n_genes_per_strain
        infection.genes[:,i] = s.gene_pool[:, rand(1:size(s.gene_pool)[2])]
    end
    
    # Add infection to host
    push!(host.liver_infections, infection)
    
    true
end


### SWITCHING EVENT ###

function get_rate_switching(t, s)
    P.switching_rate * P.n_hosts * s.n_active_infections_per_host_max
end

function do_switching!(t, s, stats)
    index = rand(CartesianIndices((P.n_hosts, s.n_active_infections_per_host_max)))
    host = s.hosts[index[1]]
    inf_index = index[2]
    
    # Advance host (rebirth or infection activation)
    advance_host!(t, s, host)
    
    # If the infection index is out of range, this is a rejected sample.
    # Otherwise we'll proceed.
    if inf_index > length(host.active_infections)
        return false
    end
    infection = host.active_infections[inf_index]
    
    # Advance expression until a non-immune gene is reached
    while true
        # Increment immunity level to currently expressed gene
        increment_immunity!(s, host, infection.genes[:, infection.expression_index])
        
        # If we're at the end, clear the infection and return
        if infection.expression_index == P.n_genes_per_strain
            delete_and_swap_with_end!(host.active_infections, inf_index)
            return true
        end
        
        # Otherwise, advance expression
        infection.expression_index += 1
        
        # If the host not immune, stop advancing
        if !is_immune(host, infection.genes[:, infection.expression_index])
            return true
        end
    end
end


### MUTATION EVENT ###

function get_rate_mutation(t, s)
    P.mutation_rate * P.n_hosts * s.n_active_infections_per_host_max * P.n_genes_per_strain * P.n_loci
end

function do_mutation!(t, s, stats)
    index = rand(CartesianIndices((P.n_hosts, s.n_active_infections_per_host_max, P.n_genes_per_strain, P.n_loci)))
    host = s.hosts[index[1]]
    inf_index = index[2]
    expression_index = index[3]
    locus = index[4]
    
    # Advance host (rebirth or infection activation)
    advance_host!(t, s, host)
    
    # If there's no active infection at the drawn index, reject this sample.
    if inf_index > length(host.active_infections)
        return false
    end
    
    infection = host.active_infections[inf_index]
    
#     println("do_mutation!($(t))")
#     println("host = $(host.id), inf = $(infection.id), locus = $(locus)")
    
    # If we ever generate too many alleles for 16-bit ints, we'll need to use bigger ones.
    @assert s.n_alleles[locus] < typemax(AlleleId)
    
    # Generate a new allele and insert it at the drawn location.
    s.n_alleles[locus] += 1
    infection.genes[locus, expression_index] = s.n_alleles[locus]
    infection.strain_id = next_strain_id!(s)
    
    true
end


### ECTOPIC RECOMBINATION EVENT ###

function get_rate_ectopic_recombination(t, s)
    P.ectopic_recombination_rate *
        P.n_hosts * s.n_active_infections_per_host_max *
        P.n_genes_per_strain * (P.n_genes_per_strain - 1) / 2.0
end

function do_ectopic_recombination!(t, s, stats)
    index = rand(CartesianIndices((P.n_hosts, s.n_active_infections_per_host_max)))
    host = s.hosts[index[1]]
    inf_index = index[2]
    
    # Advance host (rebirth or infection activation)
    advance_host!(t, s, host)
    
    # If there's no active infection at the drawn index, reject this sample.
    if inf_index > length(host.active_infections)
        return false
    end
    
    infection = host.active_infections[inf_index]
    
    gene_index_1 = rand(1:P.n_genes_per_strain)
    gene_index_2 = rand(1:P.n_genes_per_strain)
    
    gene1 = infection.genes[:, gene_index_1]
    gene2 = infection.genes[:, gene_index_2]
    
    # If the genes are the same, this is a no-op
    if gene1 == gene2
        return false
    end
    
    # Choose a breakpoint
    breakpoint = rand(1:P.n_loci)
    if breakpoint == 1
        return false
    end
    
    p_viable = p_recombination_is_viable(gene1, gene2, breakpoint)
    
#     println("do_ectopic_recombination!($(t))")
#     println("host = $(host.id), inf = $(infection.id), breakpoint = $(breakpoint), p_viable = $(p_viable)")
    
    recombined = false
    
    # Recombine to modify first gene, if viable
    if rand() < p_viable
        infection.genes[:, gene_index_1] = recombine_genes(gene1, gene2, breakpoint)
        recombined = true
    end
    
    # Recombine to modify second gene, if viable
    if rand() < p_viable
        infection.genes[:, gene_index_2] = recombine_genes(gene2, gene1, breakpoint)
        recombined = true
    end
    
    if recombined
        infection.strain_id = next_strain_id!(s)
        true
    else
        false
    end
end

"""
    Model for probability that a recombination is viable.
    
    TODO: detailed description
"""
function p_recombination_is_viable(gene1, gene2, breakpoint)
    n_diff = 0 # was "p_div"
    n_diff_before = 0 # was "child_div"
    rho = 0.8
    avg_mutation = 5.0
    for i in 1:P.n_loci
        if gene1[i] != gene2[i]
            n_diff += 1
            if i < breakpoint
                n_diff_before += 1
            end
        end
    end
    rho_power = n_diff_before * avg_mutation *
        (n_diff - n_diff_before) * avg_mutation /
        (n_diff * avg_mutation - 1.0)

    rho^rho_power
end

function recombine_genes(gene1, gene2, breakpoint)
    gene = MGene(undef)
    gene[1:(breakpoint - 1)] = gene1[1:(breakpoint - 1)]
    gene[breakpoint:end] = gene2[breakpoint:end]
    gene
end


### IMMUNITY LOSS EVENT ###

function get_rate_immunity_loss(t, s)
    P.immunity_loss_rate * P.n_hosts * s.n_immunities_per_host_max
end

function do_immunity_loss!(t, s, stats)
    index = rand(CartesianIndices((P.n_hosts, s.n_immunities_per_host_max)))
    host = s.hosts[index[1]]
    immunity_index = index[2]
    
    # Advance host (rebirth or infection activation)
    advance_host!(t, s, host)
    
    # If the immunity index is beyond this host's immunity count, reject this sample.
    if immunity_index > length(host.immunity)
        return false
    end
    
    gene = get_key_by_iteration_order(host.immunity, immunity_index)
    decrement_immunity!(host, gene)
    
#     println("do_immunity_loss($(t))")
#     println("host: $(host.id), gene: $(gene)")
    
    false
end

### MISCELLANEOUS FUNCTIONS ###

function next_host_id!(s)
    id = s.next_host_id
    s.next_host_id += 1
    id
end

function next_infection_id!(s)
    id = s.next_infection_id
    s.next_infection_id += 1
    id
end

function next_strain_id!(s)
    id = s.next_strain_id
    s.next_strain_id += 1
    id
end

function increment_immunity!(s, host, gene)
    # Get old immunity level from immunity dict
    old_level = get(host.immunity, gene, ImmunityLevel(0))
    
    # Increment immunity if the level is not at the maximum value (255 = 0xFF)
    if old_level < typemax(ImmunityLevel)
        host.immunity[gene] = old_level + 1
    end
    
    s.n_immunities_per_host_max = max(s.n_immunities_per_host_max, length(host.immunity))
end

function decrement_immunity!(host, gene)
    # Get old immunity level from immunity dict
    old_level = get(host.immunity, gene, ImmunityLevel(0))
    
    # Decrement immunity, removing it entirely if we reach 0
    if old_level == 1
        delete!(host.immunity, gene)
    else
        host.immunity[gene] -= 1
    end
end

function is_immune(host, gene)
    get(host.immunity, gene, 0) > 0
end

