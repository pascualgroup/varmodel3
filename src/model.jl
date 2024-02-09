"""
This file loads code shared by all model implementations, and then loads the
variant of the model specified by the global parameters constant P.

The `run()` function, called in run scripts, is defined by the specific variant
of the model that is loaded, e.g., inside `continuous/model.jl`.

This is an unsophisticated but straightforward way to have multiple variants of
the model share a parameters format, run script, and other bits of code.
"""

# Verify that "parameters.jl" got loaded already, and that the global parameters
# constant P is defined
@assert @isdefined Params
@assert @isdefined P
@assert typeof(P) === Params
validate(P)

include("util.jl")
include("state.jl")
include("output.jl")

const N_EVENTS = 6
const EVENTS = collect(1:N_EVENTS)
const (BITING, IMMIGRATION, SWITCHING, MUTATION, ECTOPIC_RECOMBINATION, IMMUNITY_LOSS) = EVENTS

function run()
    db = initialize_database()

    # Seed the random number generator using the provided seed,
    # or, if absent, by generating one from the OS's source of entropy.
    rng_seed = if isnothing(P.rng_seed)
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
        event_happened = do_event!(t, s, stats, event, db)

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
    s.n_immunities_per_host_max = maximum(immuneLength(host.immunity) for host in s.hosts)
    s.n_active_infections_per_host_max = maximum(length(host.active_infections) for host in s.hosts)
    s.n_liver_infections_per_host_max = maximum(length(host.liver_infections) for host in s.hosts)
end


### INITIALIZATION ###

function initialize_state()
    println("initialize_state()")

    # Initialize gene pool as an (n_loci, n_genes_initial) matrix filled with
    # allele IDs drawn uniformly randomly in 1:n_alleles_per_locus_initial.
    # and to avoid duplications, 10 times of random numbers are draw and then reduced to a set
    gene_pool = reshape(
        rand(1:P.n_alleles_per_locus_initial, P.n_loci * P.n_genes_initial*10),
        (P.n_loci, P.n_genes_initial*10)
    )
    gene_pool_set = Set()
    i = 1
    while length(gene_pool_set)<P.n_genes_initial
        push!(gene_pool_set, gene_pool[:,i])
        i+=1
    end
    gene_pool = Array{Int}(undef, P.n_loci, 0)
    for gene in gene_pool_set
        gene_pool = hcat(gene_pool, gene)
    end

    # Initialize n_hosts hosts, all born at t = 0, with lifetime drawn from a
    # distribution, and no initial infections or immunity.
    ImmuneHistoryType = if P.use_immunity_by_allele
        ImmuneHistoryByAllele
    else
        ImmuneHistoryByGene
    end
    hosts = [
        Host(
            id = id,
            t_birth = 0, t_death = draw_host_lifetime(),
            liver_infections = [], active_infections = [],
            immunity = ImmuneHistoryType(),
            n_cleared_infections = 0
        )
        for id in 1:P.n_hosts
    ]
    for host in hosts
        host.t_birth = -(rand() * host.t_death)
        host.t_death = host.t_death + host.t_birth
    end

    # Infect n_initial_infections hosts at t = 0. Genes in the infection
    # are sampled uniformly randomly from the gene pool.
    for (infection_id, host_index) in enumerate(sample(1:P.n_hosts, P.n_initial_infections, replace = false))
        infection = create_empty_infection()
        infection.id = infection_id
        infection.t_infection = 0.0
        infection.t_expression = NaN
        infection.t_last_switch = NaN
        infection.duration = NaN
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
        n_active_infections_per_host_max = 0,
        n_liver_infections_per_host_max = 0,
        n_cleared_infections = 0,
        durations = []
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
        t_expression = NaN,
        t_last_switch = NaN,
        duration = NaN,
        strain_id = StrainId(0),
        genes = fill(AlleleId(0), (P.n_loci, P.n_genes_per_strain)),
        expression_index = ExpressionIndex(0),
        expression_index_locus = ExpressionIndexLocus(0)
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

function do_event!(t, s, stats, event, db)
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

    stats.n_bites += 1

    # Uniformly randomly sample infecting host (source) and host being infected
    # (destination).
    src_host = rand(s.hosts)
    dst_host = rand(s.hosts)

    # Advance host (rebirth or infection activation)
    advance_host!(t, s, src_host, stats)
    advance_host!(t, s, dst_host, stats)

    # The source host must be infected in order to transmit.
    src_active_count = length(src_host.active_infections)
    if src_active_count == 0
        return false
    end
    stats.n_infected_bites += 1

    # The destination host must have space available in the liver stage.
    dst_available_count = if isnothing(P.n_infections_liver_max)
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

    # First choose the active strains from the host that will be transmitted to mosquito
    # This is determined by the transmissibility
    choose_transmit = rand(Float64, src_active_count)
    transmitted_strains = src_host.active_infections[choose_transmit.<p_transmit]
    #println("t = $(t): originalSize $(src_active_count), newSize $(length(transmitted_strains))")

    # The number of transmissions is bounded by the number of source infections
    # and the number of available slots in the destination.
    n_transmissions_max = min(length(transmitted_strains), dst_available_count)
    transmitted = false
    for i in 1:n_transmissions_max
        
        stats.n_transmissions += 1
        transmitted = true

        # Randomly sample two source infections within transmitted_strains to recombine
        src_inf_1 = rand(transmitted_strains)
        src_inf_2 = rand(transmitted_strains)

        # Get a new infection struct, or recycle an old infection
        # to prevent excess memory allocation.
        dst_inf = recycle_or_create_infection(s)
        dst_inf.id = next_infection_id!(s)
        dst_inf.t_infection = t
        dst_inf.t_expression = NaN
        dst_inf.expression_index = 0
        dst_inf.expression_index_locus = 0
        dst_inf.duration = NaN

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
        # update population wide host liver max
        s.n_liver_infections_per_host_max = max(s.n_liver_infections_per_host_max, length(dst_host.liver_infections))

    end

    if transmitted
        stats.n_transmitting_bites += 1
        true
    else
        false
    end
end

function advance_host!(t, s, host, stats)
    if t > host.t_death
        # If the host is past its death time
        do_rebirth!(t, s, host)
    else
        i = 1
        while i <= length(host.liver_infections)
            infection = host.liver_infections[i]
            if infection.t_infection + P.t_liver_stage < t
                # println("t = $(t): activating host $(host.id), inf $(infection.id)")

                # If the infection is past the liver stage, remove it from the liver.
                delete_and_swap_with_end!(host.liver_infections, i)
                # If there's room, move it into the active infections array.
                # Otherwise, just put it into the recycle bin.
                if isnothing(P.n_infections_active_max) || length(host.active_infections) < P.n_infections_active_max
                    infection.expression_index = 1
                    if P.whole_gene_immune
                        infection.expression_index_locus = P.n_loci
                    else
                        infection.expression_index_locus = 1
                    end
                    push!(host.active_infections, infection)
                    infection.t_expression = t
                    infection.t_last_switch = t

                    advance_immuned_genes!(t,s,host,length(host.active_infections), stats)

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
    host.n_cleared_infections = 0
    empty!(host.liver_infections)
    empty!(host.active_infections)
    empty!(host.immunity)
end


### IMMIGRATION EVENT ###

function get_rate_immigration(t, s)
    P.immigration_rate_fraction * get_rate_biting(t, s)
end

function do_immigration!(t, s, stats)

    # Sample a random host and advance it (rebirth or infection activation)
    host = rand(s.hosts)
    advance_host!(t, s, host, stats)

    # If host doesn't have an available infection slot, reject this sample.
    if isnothing(P.n_infections_liver_max)
        if length(host.liver_infections) == P.n_infections_liver_max
            return false
        end
    end

    # Construct infection by sampling from gene pool
    infection = recycle_or_create_infection(s)
    infection.id = next_infection_id!(s)
    infection.t_infection = t
    infection.t_expression = NaN
    infection.duration = NaN
    infection.strain_id = next_strain_id!(s)
    infection.expression_index = 0
    infection.expression_index_locus = 0
    for i in 1:P.n_genes_per_strain
        infection.genes[:,i] = s.gene_pool[:, rand(1:size(s.gene_pool)[2])]
    end

    # Add infection to host
    push!(host.liver_infections, infection)

    s.n_liver_infections_per_host_max = max(s.n_liver_infections_per_host_max, length(host.liver_infections))

    true
end


### SWITCHING EVENT ###

function get_rate_switching(t, s)
    if P.use_immunity_by_allele && !P.whole_gene_immune
        #switching rate set by total number of alleles
        (P.switching_rate * P.n_loci) * P.n_hosts * s.n_active_infections_per_host_max
    else
        P.switching_rate * P.n_hosts * s.n_active_infections_per_host_max
    end
end

function do_switching!(t, s, stats)
    # change to total number of infections instead of active infections alone
    index = rand(CartesianIndices((P.n_hosts, s.n_active_infections_per_host_max)))
    host = s.hosts[index[1]]
 
    # Advance host (rebirth or infection activation)
    advance_host!(t, s, host, stats)
    inf_index = index[2]

    # If the infection index is out of range, this is a rejected sample.
    # Otherwise we'll proceed.
    if inf_index > length(host.active_infections)
        return false
    end
    infection = host.active_infections[inf_index]

    """
    Increment immunity level to currently expressed gene.
    For the partial allele model, expression advance by alleles, but
    Immunity only gains after the full gene finishes expression
    """
    if infection.expression_index_locus == P.n_loci
        increment_immunity!(s, host, infection.genes[:, infection.expression_index])
    end
    
    stats.n_switch_not_immune += 1
    stats.t_switch_sum += t - infection.t_last_switch
    
    # If we're at the end, clear the infection and return.
    if infection.expression_index == P.n_genes_per_strain && infection.expression_index_locus == P.n_loci
        if s.n_cleared_infections % P.sample_duration == 0
            # Calculate and write the infection duration.
            add_infectionDuration!(t, s, host, inf_index)
        end
        s.n_cleared_infections += 1
        host.n_cleared_infections += 1
        delete_and_swap_with_end!(host.active_infections, inf_index)
    else
        # Otherwise, advance gene and/or locus expression(s).
        if P.use_immunity_by_allele  && !P.whole_gene_immune
            if infection.expression_index_locus == P.n_loci
                infection.expression_index += 1
                infection.expression_index_locus = 1
            else
                infection.expression_index_locus += 1
            end
        else
            infection.expression_index += 1
            infection.expression_index_locus = P.n_loci
        end
        infection.t_last_switch = t
    end
    
    """
    after gaining immunity to the expressed gene, loop through the other infections
    in the host to see if any one needs advancing
    """
    i = 1
    while i <= length(host.active_infections)
        if advance_immuned_genes!(t, s, host, i, stats)
            # if there is no end of expression and reordering of infections, then index plus 1
            i += 1
        end
    end
 
    return true
end

# This function moves the expression index of an infection to its first non-immune allele/gene.
function advance_immuned_genes!(t, s, host, i, stats)

    # Advance expression until a non-immune gene or allele is reached.
    infection = host.active_infections[i]
    # If the host not immune, stop advancing.
    to_advance = is_immune(host.immunity, infection.genes[:, infection.expression_index],infection.expression_index_locus)
 
    while to_advance
        # Record an instantaneous switching event
        stats.n_switch_immune += 1
        
        # Increment immunity level to currently expressed gene or allele.
        if infection.expression_index_locus == P.n_loci
            increment_immunity!(s, host, infection.genes[:, infection.expression_index])
        end
 
        # If we're at the end, clear the infection and return.
        if infection.expression_index == P.n_genes_per_strain && infection.expression_index_locus == P.n_loci
            if s.n_cleared_infections % P.sample_duration == 0 
                # Calculate and write the infection duration.
                add_infectionDuration!(t, s, host, i)
            end
            s.n_cleared_infections += 1
            host.n_cleared_infections += 1
            delete_and_swap_with_end!(host.active_infections, i)
            #here if deleting an infection, and the indexing changed, then tell the calling function 
            #it doesn't advancing its index
            return false
        else
            # Otherwise, advance gene and/or locus expression(s).
            if P.use_immunity_by_allele  && !P.whole_gene_immune
                if infection.expression_index_locus == P.n_loci
                    infection.expression_index += 1
                    infection.expression_index_locus = 1
                else
                    infection.expression_index_locus += 1
                end
            else
                infection.expression_index += 1
                infection.expression_index_locus = P.n_loci
            end 
        end
        # If the host not immune, stop advancing.
        to_advance = is_immune(host.immunity, infection.genes[:, infection.expression_index],infection.expression_index_locus)
     end
     true
 end
 

### MUTATION EVENT ###
#update mutation and recombinatin rates towards all infections
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
    advance_host!(t, s, host, stats)

    # If there's no active infection at the drawn index, reject this sample.
    if inf_index > length(host.active_infections)
        return false
    end

    infection = host.active_infections[inf_index]

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
    advance_host!(t, s, host, stats)

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

    breakpoint, p_functional = if P.ectopic_recombination_generates_new_alleles
        # Choose a breakpoint
        breakpoint = P.n_loci * rand()
        p_functional = p_recombination_is_functional_real(gene1, gene2, breakpoint)

        (Int(ceil(breakpoint)), p_functional)
    else
        # Choose a breakpoint
        breakpoint = rand(1:P.n_loci)
        if breakpoint == 1
            return false
        end

        p_functional = p_recombination_is_functional_integer(gene1, gene2, breakpoint)

        (breakpoint, p_functional)
    end

    is_conversion = rand() < P.p_ectopic_recombination_is_conversion

    recombined = false

    create_new_allele = P.ectopic_recombination_generates_new_alleles &&
        rand() < P.p_ectopic_recombination_generates_new_allele

    # Recombine to modify first gene, if functional
    if !is_conversion && rand() < p_functional
        infection.genes[:, gene_index_1] = if create_new_allele
            recombine_genes_new_allele(s, gene1, gene2, breakpoint)
        else
            recombine_genes(gene1, gene2, breakpoint)
        end
        recombined = true
    end

    # Recombine to modify second gene, if functional
    if rand() < p_functional
        infection.genes[:, gene_index_2] = if create_new_allele
            recombine_genes_new_allele(s, gene1, gene2, breakpoint)
        else
            recombine_genes(gene2, gene1, breakpoint)
        end
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
    Probability that a recombination results in a functional gene.

    Version for real-valued breakpoint, used when
    `ectopic_recombination_generates_new_alleles == true`.
"""
function p_recombination_is_functional_real(gene1, gene2, breakpoint::Float64)
    n_diff = 0.0 # was "p_div"
    n_diff_before = 0.0 # was "child_div"
    rho = P.rho_recombination_tolerance
    mean_n_mutations = P.mean_n_mutations_per_epitope
    for i in 1:P.n_loci
        if gene1[i] != gene2[i]
            n_diff += 1
            if i - 1 < breakpoint
                if breakpoint - i > 0
                    n_diff_before += 1
                else
                    n_diff_before += breakpoint - (i - 1)
                end
            end
        end
    end
    rho_power = n_diff_before * mean_n_mutations *
        (n_diff - n_diff_before) * mean_n_mutations /
        (n_diff * mean_n_mutations - 1.0)

    rho^rho_power
end

"""
    Model for probability that a recombination results in a functional gene.

    Version for integer-valued breakpoint, used when
    `ectopic_recombination_generates_new_alleles == false`.
"""
function p_recombination_is_functional_integer(gene1, gene2, breakpoint::Int)
    n_diff = 0 # was "p_div"
    n_diff_before = 0 # was "child_div"
    rho = P.rho_recombination_tolerance
    mean_n_mutations = P.mean_n_mutations_per_epitope
    for i in 1:P.n_loci
        if gene1[i] != gene2[i]
            n_diff += 1
            if i < breakpoint
                n_diff_before += 1
            end
        end
    end
    rho_power = n_diff_before * mean_n_mutations *
        (n_diff - n_diff_before) * mean_n_mutations /
        (n_diff * mean_n_mutations - 1.0)

    rho^rho_power
end

function recombine_genes(gene1, gene2, breakpoint)
    gene = MGene(undef)
    gene[1:(breakpoint - 1)] = gene1[1:(breakpoint - 1)]
    gene[breakpoint:end] = gene2[breakpoint:end]
    gene
end

function recombine_genes_new_allele(s, gene1, gene2, breakpoint)
    gene = MGene(undef)
    gene[1:(breakpoint - 1)] = gene1[1:(breakpoint - 1)]
    if gene1[breakpoint] != gene2[breakpoint]
        # If we ever generate too many alleles for 16-bit ints, we'll need to use bigger ones.
        @assert s.n_alleles[breakpoint] < typemax(AlleleId)

        s.n_alleles[breakpoint] += 1
        gene[breakpoint] = s.n_alleles[breakpoint]
    else
        gene[breakpoint] = gene2[breakpoint]
    end
    if P.n_loci > breakpoint
        gene[(breakpoint + 1):end] = gene2[(breakpoint + 1):end]
    end
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
    advance_host!(t, s, host, stats)

    # If the immunity index is beyond this host's immunity count, reject this sample.
    if immunity_index >  immuneLength(host.immunity)
        return false
    end

    decrement_immunity_at_sampled_index!(host.immunity, immunity_index)

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


### IMMUNITY FUNCTIONS ###

function empty!(ih::ImmuneHistoryByGene)
    empty!(ih.d)
end

function empty!(ih::ImmuneHistoryByAllele)
    for d in ih.vd
        empty!(d)
    end
end


function increment_immunity!(s, host, gene)
    increment_immunity!(host.immunity, gene)
    s.n_immunities_per_host_max = max(s.n_immunities_per_host_max, immuneLength(host.immunity))
end

function increment_immunity!(ih::ImmuneHistoryByGene, gene)
    # Get old immunity level from immunity dict
    old_level = get(ih.d, gene, ImmunityLevel(0))

    # Increment immunity if the level is not at the maximum value (255 = 0xFF)
    if old_level < typemax(ImmunityLevel)
        ih.d[gene] = old_level + 1
    end
end

function increment_immunity!(ih::ImmuneHistoryByAllele, gene)
    # Increment immunity at each locus
    for (locus, allele_id) in enumerate(gene)
        old_level = get(ih.vd[locus], allele_id, ImmunityLevel(0))
        if old_level < typemax(ImmunityLevel)
            ih.vd[locus][allele_id] = old_level + 1
        end
    end
end


function immuneLength(ih::ImmuneHistoryByGene)
    length(ih.d)
end

function immuneLength(ih::ImmuneHistoryByAllele)
    sum(length(ih.vd[locus]) for locus in 1:length(ih.vd))
end


function decrement_immunity_at_sampled_index!(ih::ImmuneHistoryByGene, index)
    gene = get_key_by_iteration_order(ih.d, index)
    decrement_immunity!(ih, gene)
end

function decrement_immunity!(ih::ImmuneHistoryByGene, gene)
    # Get old immunity level from immunity dict
    old_level = ih.d[gene]

    # Decrement immunity, removing it entirely if we reach 0
    if old_level == 1
        delete!(ih.d, gene)
    else
        ih.d[gene] -= 1
    end
end

function decrement_immunity_at_sampled_index!(ih::ImmuneHistoryByAllele, index)
    cur_index = index
    for locus in 1:P.n_loci
        if cur_index <= length(ih.vd[locus])
            allele_id = get_key_by_iteration_order(ih.vd[locus], cur_index)
            decrement_immunity!(ih, Locus(locus), allele_id)
            break
        else
            cur_index -= length(ih.vd[locus])
        end
    end
end

function decrement_immunity!(ih::ImmuneHistoryByAllele, locus::Locus, allele_id::AlleleId)
    # Get old immunity level from immunity dict.
    old_level = ih.vd[locus][allele_id]

    # Decrement immunity, removing it entirely if we reach 0.
    if old_level == 1
        delete!(ih.vd[locus], allele_id)
    else
        ih.vd[locus][allele_id] -= 1
    end
end

function is_immune(ih::ImmuneHistoryByGene, gene, loc)
    get(ih.d, gene, 0) > 0
end

function is_immune(ih::ImmuneHistoryByAllele, gene, loc)
    if P.whole_gene_immune
        for (locus, allele_id) in enumerate(gene)
            if get(ih.vd[locus], allele_id, 0) == 0
                return false
            end
        end
    else
        if get(ih.vd[loc], gene[loc], 0) == 0
            return false
        end
    end
    true
end

# push the current infection's duration calculations into durations vector
function add_infectionDuration!(t, s, host, i)
    get_duration!(host.active_infections, i, t)
    newInfDur = infectionDuration(
        id=host.active_infections[i].id,
        host_id=host.id,
        n_cleared_infections=host.n_cleared_infections,
        n_immuned_alleles=immuneLength(host.immunity),
        t_infection = host.active_infections[i].t_infection,
        t_expression = host.active_infections[i].t_expression,
        duration = host.active_infections[i].duration
    )
    push!(s.durations, newInfDur)
end
