include("util.jl")
include("state.jl")
include("output.jl")

const N_EVENTS = 5
const EVENTS = collect(1:N_EVENTS)
const (BITING, IMMIGRATION, SWITCHING, MUTATION, ECTOPIC_RECOMBINATION) = EVENTS

function run_exact()
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
    db = initialize_database()
    write_output(0, s, db)
    
    # Initialize event rates
    rates = fill(0.0, N_EVENTS)
    total_rate, dt_dist = update_rates!(t, s, rates)
    
    # Loop events until end of simulation
    while t < P.t_end
        # Draw next time with rate equal to the sum of all event rates
        dt = rand(dt_dist)
        event = direct_sample_linear_scan(rates, total_rate)
        
        # At each integer time, write output/state verification (if necessary)
        # and update rates to ensure that biting rate is kept current.
        # Loop required in case the simulation jumps past two integer times.
        while t_next_integer < t + dt
            t_next = Float64(t_next_integer)
            write_output(t_next, s, db)
            if P.verification_period != nothing && t_next_integer % P.verification_period == 0
                verify(t_next, s)
            end
            total_rate, dt_dist = update_rates!(t, s, rates)
            t_next_integer += 1
        end
        
        # Update time
        t += dt
        
        # Execute event.
        # Many events are no-ops due to rejection sampling.
        # If an event actually happened, update event rates.
        event_happened = do_event!(t, s, event)
        if event_happened
            total_rate, dt_dist = update_rates!(t, s, rates)
        end
    end
    
    println("elapsed time (s): $(Dates.value(now() - last_summary_datetime) / 1000.0)")
end


### INITIALIZATION ###

function initialize_state()
    gene_pool = reshape(
        rand(1:P.n_alleles_per_locus_initial, P.n_loci * P.n_genes_initial),
        (P.n_loci, P.n_genes_initial)
    )
    
    hosts = [initialize_host(id) for id in 1:P.n_hosts]
    
    # Infect n_initial_infections hosts
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

function initialize_host(id)
    lifetime = draw_host_lifetime()
    t_birth = -rand() * lifetime
    t_death = t_birth + lifetime
    
    liver_infections = []
    sizehint!(liver_infections, P.n_infections_liver_max)
    
    active_infections = []
    sizehint!(active_infections, P.n_infections_active_max)
    
    Host(
        id = id,
        t_birth = t_birth,
        t_death = t_death,
        liver_infections = liver_infections,
        active_infections = active_infections,
        immunity = Dict()
    )
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

function update_rates!(t, s, rates)
    for event in EVENTS
        rates[event] = get_rate(t, s, event)
    end
    total_rate = sum(rates)
    (total_rate, Exponential(1.0 / total_rate))
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

function do_event!(t, s, event)
    if event == BITING
        do_event_biting!(t, s)
    elseif event == IMMIGRATION
        do_event_immigration!(t, s)
    elseif event == SWITCHING
        do_event_switching!(t, s)
    elseif event == MUTATION
        do_event_mutation!(t, s)
    elseif event == ECTOPIC_RECOMBINATION
        do_event_ectopic_recombination!(t, s)
    end
end


### BITING EVENT ###

function get_rate_biting(t, s)
    1.0
end

function do_event_biting!(t, s)
#     println("do_event_biting!()")
    true
end


### IMMIGRATION EVENT ###

function get_rate_immigration(t, s)
    1.0
end

function do_event_immigration!(t, s)
#     println("do_event_immigration!()")
    true
end


### SWITCHING EVENT ###

function get_rate_switching(t, s)
    1.0
end

function do_event_switching!(t, s)
#     println("do_event_switching!()")
    true
end


### MUTATION EVENT ###

function get_rate_mutation(t, s)
    1.0
end

function do_event_mutation!(t, s)
#     println("do_event_mutation!()")
    true
end


### ECTOPIC RECOMBINATION EVENT ###

function get_rate_ectopic_recombination(t, s)
    1.0
end

function do_event_ectopic_recombination!(t, s)
#     println("do_event_ectopic_recombination!()")
    true
end


### MISCELLANEOUS FUNCTIONS ###

function draw_host_lifetime()
    min(rand(Exponential(P.mean_host_lifetime)), P.max_host_lifetime)
end


### STATE VERIFICATION ###

function verify(t, s)
    println("verify($(t), s)")
end
