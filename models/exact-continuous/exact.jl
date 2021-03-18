include("util.jl")
include("state.jl")
include("output.jl")

const N_EVENTS = 5
const EVENTS = collect(1:N_EVENTS)
const (BITING, IMMIGRATION, SWITCHING, MUTATION, ECTOPIC_RECOMBINATION) = EVENTS

function run_exact()
    db = initialize_database()
    stats = SummaryStats()
    
    # Seed the random number generator using the provided seed,
    # or, if absent, by generating one from the OS's source of entropy.
    rng_seed = if P.rng_seed === missing
        rand(RandomDevice(), UInt64)
    else
        P.rng_seed
    end
    Random.seed!(rng_seed)
    execute(db.meta, ("rng_seed", rng_seed))
    
    # Start recording elapsed time
    start_datetime = now()
    last_summary_datetime = start_datetime
    
    # Used to decide when to do output and verification
    t_next_integer = 1
    
    # Initialize state
    t = 0.0
    s = initialize_state()
    verify(t, s)
    
    # Run initial output
    write_output(0, s, db)
    
    # Initialize event rates
    rates = fill(0.0, N_EVENTS)
    total_rate = update_rates!(t, s, rates)
    
    # Loop events until end of simulation
    while total_rate > 0.0 && t < P.t_end
        # Draw next time with rate equal to the sum of all event rates
        dt = rand(Exponential(1.0 / total_rate))
        @assert dt > 0.0 && !isinf(dt)
        event = direct_sample_linear_scan(rates, total_rate)
        
        # At each integer time, write output/state verification (if necessary).
        # Loop required in case the simulation jumps past two integer times.
        while t_next_integer < t + dt
            t_next = Float64(t_next_integer)
            write_output(t_next, s, db)
            if P.verification_period != nothing && t_next_integer % P.verification_period == 0
                verify(t_next, s)
            end
            
            t_next_integer += 1
        end
        
        # Update time
        t += dt
        
        # Execute event.
        # Many events are no-ops due to rejection sampling.
        # If an event actually happened, update event rates.
        event_happened = do_event!(t, s, stats, event)
        if event_happened
            total_rate = update_rates!(t, s, rates)
        end
    end
    
    elapsed_time = Dates.value(now() - last_summary_datetime) / 1000.0
    println("elapsed time (s): $(elapsed_time)")
    execute(db.meta, ("elapsed_time", elapsed_time))
    
    went_extinct = total_rate == 0.0
    println("went extinct? $(went_extinct)")
    execute(db.meta, ("went_extinct", Int64(went_extinct)))
end


### INITIALIZATION ###

function initialize_state()
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
        old_infections = []
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
    if !empty(s.old_infections)
        pop!(s.old_infections)
    else
        create_empty_infection()
    end
end

function update_rates!(t, s, rates)
    for event in EVENTS
        rates[event] = get_rate(t, s, event)
    end
    sum(rates)
end


### EVENT DEMUX ###

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
    end
end

function do_event!(t, s, stats, event)
    if event == BITING
        do_event_biting!(t, s, stats)
    elseif event == IMMIGRATION
        do_event_immigration!(t, s, stats)
    elseif event == SWITCHING
        do_event_switching!(t, s, stats)
    elseif event == MUTATION
        do_event_mutation!(t, s, stats)
    elseif event == ECTOPIC_RECOMBINATION
        do_event_ectopic_recombination!(t, s, stats)
    end
end


### BITING EVENT ###

function get_rate_biting(t, s)
    biting_rate = P.biting_rate[1 + Int(floor(t)) % P.t_year]
    biting_rate * P.n_hosts
end

function do_event_biting!(t, s, stats)
#     println("do_event_biting!($(t), s)")
    
    stats.n_bites += 1
    
    # Uniformly randomly sample infecting host (source) and host being infected
    # (destination).
    src_host = rand(s.hosts)
    dst_host = rand(s.hosts)
    
    # The source host must be infected in order to transmit.
    src_active_count = length(src_host.active_infections)
    if src_active_count == 0
        return false
    end
    stats.n_infected_bites += 1
    
    # The destination host must have space available in the liver stage.
    dst_available_count = P.n_infections_liver_max - length(src_host.liver_infections)
    if dst_available_count == 0
        return false
    end
    stats.n_infected_bites_with_space += 1
    
    # Compute probability of each transmission
    p_transmit = if p.coinfection_reduces_transmission
        P.transmissibility
    else
        P.transmissibility / src_inf_count
    end
    
    # The number of transmissions is bounded by the number of source infections
    # and the number of available slots in the destination.
    n_transmissions_max = min(src_active_count, dst_available_count)
    transmitted = false
    for i in 1:n_transmissions_max
        if rand() < P.transmissibility
            stats.n_transmissions += 1
            transmitted = true
            
            # Randomly sample two source infections to recombine
            src_inf_1 = rand(src_host.active_infections)
            src_inf_2 = rand(src_host.active_infections)
        
            # Get a new infection struct, or recycle an old infection
            # to prevent excess memory allocation.
            dst_inf = recycle_or_create_infection(s)
            dst_inf.t_infection = t
            dst_inf.expression_index = 0
            
            # Construct strain for new infection
            if inf1.strain_id == inf2.strain_id
                # If both infections have the same strain, then the new infection
                # is given infection 1's genes with expression order shuffled.
                dst_inf.strain_id = inf1.strain_id
                shuffle_columns_to!(dst_inf.genes, src_inf_1.genes)
            else
                # Otherwise, the new infection is given a new strain constructed by
                # taking a random sample of the genes in the two source infections.
                dst_inf.strain_id = next_strain_id(s)
                sample_columns_from_two_matrices_to!(dst_inf.genes, src_inf_1.genes, src_inf_2.genes)
            end
        
            # Add this infection to the destination host
            push!(dst_host.liver_infections, dst_inf)
        end
    end
    if transmitted
        stats.n_transmitted_bites += 1
    end
    
    true
end


### IMMIGRATION EVENT ###

function get_rate_immigration(t, s)
    0.0
end

function do_event_immigration!(t, s)
#     println("do_event_immigration!()")
    true
end


### SWITCHING EVENT ###

function get_rate_switching(t, s)
    0.0
end

function do_event_switching!(t, s)
#     println("do_event_switching!()")
    true
end


### MUTATION EVENT ###

function get_rate_mutation(t, s)
    0.0
end

function do_event_mutation!(t, s)
#     println("do_event_mutation!()")
    true
end


### ECTOPIC RECOMBINATION EVENT ###

function get_rate_ectopic_recombination(t, s)
    0.0
end

function do_event_ectopic_recombination!(t, s)
#     println("do_event_ectopic_recombination!()")
    true
end


### MISCELLANEOUS FUNCTIONS ###


### STATE VERIFICATION ###

function verify(t, s)
    println("verify($(t), s)")
end
