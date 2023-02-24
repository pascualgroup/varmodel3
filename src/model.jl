"""
This file loads code shared by all model implementations, and then loads the
variant of the model specified by the global parameters constant P.

The `run()` function, called in run scripts, is defined by the specific variant
of the model that is loaded, e.g., inside `continuous/model.jl`.

This is an unsophisticated but straightforward way to have multiple variants of
the model share a parameters format, run script, and other bits of code.
"""

# Verify that "parameters.jl" got loaded already, and that the global parameters
# constant P is defined.
@assert @isdefined Params
@assert @isdefined P
@assert typeof(P) === Params
validate(P)

# Parameters that need (vector-of-vector)-to-matrix conversion
if P.snp_linkage_disequilibrium
    const P_snp_pairwise_ld = reduce(hcat, P.snp_pairwise_ld)
end

include("util.jl")
include("state.jl")
include("output.jl")

if P.generalized_immunity
    const N_EVENTS = 8
    const EVENTS = collect(1:N_EVENTS)
    const (BITING, IMMIGRATION, BACKGROUND_CLEARANCE, SWITCHING, MUTATION, ECTOPIC_RECOMBINATION, IMMUNITY_LOSS, GENERALIZED_IMMUNITY_LOSS) = EVENTS
else
    const N_EVENTS = 7
    const EVENTS = collect(1:N_EVENTS)
    const (BITING, IMMIGRATION, BACKGROUND_CLEARANCE, SWITCHING, MUTATION, ECTOPIC_RECOMBINATION, IMMUNITY_LOSS) = EVENTS
end

const USE_BITING_RATE_MULTIPLIER_BY_YEAR = P.biting_rate_multiplier_by_year !== nothing

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

    # Start recording elapsed time.
    start_datetime = now()

    # Used to decide when to do output and verification.
    t_next_integer = 1

    # Initialize state.
    t = 0.0
    s = initialize_state()
    verify(t, s)

    # Run initial output.
    stats = SummaryStats()
    write_output!(db, 0, s, stats)

    if P.n_snps_per_strain > 0
        write_initial_snp_allele_frequencies(db, s)
    end

    # Initialize event rates.
    rates = [get_rate(t, s, event) for event in EVENTS]
    total_rate = sum(rates)

    # Loop events until end of simulation.
    while total_rate > 0.0 && t < P.t_end
        # Draw next time with rate equal to the sum of all event rates.
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

            if P.migrants_match_local_prevalence
                if t_next_integer % P.migration_rate_update_period == 0
                    recompute_infected_ratio!(s)
                end
            end

            t_next_integer += 1
        end

        # Draw the event, update time, and execute event.
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

    total_num_genes_generated_mut = s.next_gene_id_mut - 1
    println("total number of new genes generated out of mutation? $(total_num_genes_generated_mut)")
    execute(db.meta, ("total_num_genes_generated_mut", Int64(total_num_genes_generated_mut)))

    total_num_genes_generated_recomb = s.next_gene_id_recomb - 1
    println("total number of new genes generated out of recombination? $(total_num_genes_generated_recomb)")
    execute(db.meta, ("total_num_genes_generated_recomb", Int64(total_num_genes_generated_recomb)))
end

function recompute_rejection_upper_bounds!(s)
    s.n_immunities_per_host_max = maximum(immunity_count(host.immunity) for host in s.hosts)
    s.n_active_infections_per_host_max = maximum(length(host.active_infections) for host in s.hosts)
    s.n_liver_infections_per_host_max = maximum(length(host.liver_infections) for host in s.hosts)
    println("recompute_rejection_upper_bounds!: s.n_active_infections_per_host_max $(s.n_active_infections_per_host_max)")
    println("recompute_rejection_upper_bounds!: s.n_liver_infections_per_host_max $(s.n_liver_infections_per_host_max)")
end

function recompute_infected_ratio!(s)
    s.infected_ratio = if s.n_bites_for_migration_rate == 0
        0.0
    else
        s.n_transmitting_bites_for_migration_rate / s.n_bites_for_migration_rate
    end
    s.n_transmitting_bites_for_migration_rate = 0
    s.n_bites_for_migration_rate = 0
end


### INITIALIZATION ###

#function initialize_state(a)
function initialize_state()
    println("initialize_state()")

    # Initialize gene pool as an (n_loci, n_genes_initial) matrix whose columns
    # are unique randomly generated genes.
    gene_pool_set = Set()
    while length(gene_pool_set) < P.n_genes_initial
        push!(gene_pool_set, rand(1:P.n_alleles_per_locus_initial, P.n_loci))
    end
    gene_pool = zeros(AlleleId, P.n_loci, P.n_genes_initial)
    for (i, gene) in enumerate(gene_pool_set)
        gene_pool[:,i] = gene
    end

    # Initialize n_hosts hosts, all born at t = 0, with lifetime drawn from a
    # distribution, and no initial infections or immunity.
    hosts = [
        Host(
            id = id,
            t_birth = 0, t_death = draw_host_lifetime(),
            liver_infections = [], active_infections = [], active_infections_detectable = [],
            immunity = ImmuneHistory(), generalized_immunity = 0,
            n_cleared_infections = 0
        )
        for id in 1:P.n_hosts
    ]
    for host in hosts
        host.t_birth = -(rand() * host.t_death)
        host.t_death = host.t_death + host.t_birth
    end

    # Define the initial allele frequencies at each SNP.
    snp_allele_freq = fill(0.5, P.n_snps_per_strain)
    if P.n_snps_per_strain > 0
        if P.distinct_initial_snp_allele_frequencies
            if !P.snp_linkage_disequilibrium
                # As unlinked SNPs are independent, their initial allele frequencies
                # are independently defined.
                for snp in 1:P.n_snps_per_strain
                    snp_allele_freq[snp] = rand(Uniform(P.initial_snp_allele_frequency[1], P.initial_snp_allele_frequency[2]))
                    snp_allele_freq[snp] = rand([snp_allele_freq[snp], 1 - snp_allele_freq[snp]])
                end
            else
                unlinked_snps = collect(1:P.n_snps_per_strain)
                i = 1
                while i in unlinked_snps
                    # Define the linked SNPs:
                    linked_snps = find_linked_snps(i)

                    # As linked SNPs are related, their initial allele frequencies
                    # are co-defined.
                    if size(linked_snps)[1] > 1
                        unlinked_snps = setdiff(unlinked_snps, linked_snps)
                        snp_allele_freq[linked_snps[1]] = rand(Uniform(P.initial_snp_allele_frequency[1], P.initial_snp_allele_frequency[2]))
                        snp_allele_freq[linked_snps[1]] = rand([snp_allele_freq[linked_snps[1]], 1 - snp_allele_freq[linked_snps[1]]])
                        for linked_snp in linked_snps[2:size(linked_snps)[1]]
                            snp_allele_freq[linked_snp] = snp_allele_freq[linked_snps[1]]
                            snp_allele_freq[linked_snp] = rand([snp_allele_freq[linked_snp], 1 - snp_allele_freq[linked_snp]])
                        end
                    end
                    i += 1
                    while !(i in unlinked_snps) && i <= P.n_snps_per_strain
                        i += 1
                    end
                end

                # As unlinked SNPs are independent, their initial allele frequencies
                # are independently defined.
                for unlinked_snp in unlinked_snps
                    snp_allele_freq[unlinked_snp] = rand(Uniform(P.initial_snp_allele_frequency[1], P.initial_snp_allele_frequency[2]))
                    snp_allele_freq[unlinked_snp] = rand([snp_allele_freq[unlinked_snp], 1 - snp_allele_freq[unlinked_snp]])
                end
            end
        end
    end

    # Infect n_initial_infections hosts at t = 0. Genes in the infection
    # are sampled uniformly randomly from the gene pool.
    for (infection_id, host_index) in enumerate(sample(1:P.n_hosts, P.n_initial_infections, replace = false))
        infection = create_empty_infection()
        infection.id = infection_id
        infection.t_infection = 0.0
        infection.t_expression = NaN
        infection.duration = NaN
        if P.drug_treatment
            infection.p_symptoms = P.pathogenicity
        end
        infection.strain_id = infection_id
        infection.genes[:,:] = reshape(
            gene_pool[:, rand(1:P.n_genes_initial, P.n_genes_per_strain)],
            (P.n_loci, P.n_genes_per_strain)
        )
        if P.n_snps_per_strain > 0
            for snp in 1:P.n_snps_per_strain
                infection.snps[snp] = sample([1, 2],
                Weights([snp_allele_freq[snp], 1 - snp_allele_freq[snp]]))
            end
            if P.drug_treatment && P.resistant_snp && infection.snps[1] == 2
                infection.p_transmit = infection.p_transmit * P.resistant_cost
                infection.p_detect = infection.p_detect * P.resistant_cost
                infection.p_symptoms = infection.p_symptoms * P.resistant_cost
            end
        end
        push!(hosts[host_index].liver_infections, infection)
    end

    State(
        n_alleles = fill(P.n_alleles_per_locus_initial, P.n_loci),
        gene_pool = gene_pool,
        next_host_id = P.n_hosts + 1,
        next_strain_id = P.n_initial_infections + 1,
        next_infection_id = P.n_initial_infections + 1,
        next_gene_id_mut = 1,
        next_gene_id_recomb = 1,
        hosts = hosts,
        old_infections = [],
        n_immunities_per_host_max = 0,
        n_active_infections_per_host_max = 0,
        n_liver_infections_per_host_max = 0,
        n_cleared_infections = 0,
        durations = [],
        n_transmitting_bites_for_migration_rate = 0,
        n_bites_for_migration_rate = 0,
        infected_ratio = 1.0,
        initial_snp_allele_frequencies = snp_allele_freq
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
        duration = NaN,
        p_detect = P.detectability,
        p_transmit = P.transmissibility,
        p_symptoms = NaN,
        strain_id = StrainId(0),
        genes = fill(AlleleId(0), (P.n_loci, P.n_genes_per_strain)),
        snps = fill(SnpId(0), P.n_snps_per_strain),
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
    elseif event == BACKGROUND_CLEARANCE
        get_rate_background_clearance(t, s)
    elseif event == SWITCHING
        get_rate_switching(t, s)
    elseif event == MUTATION
        get_rate_mutation(t, s)
    elseif event == ECTOPIC_RECOMBINATION
        get_rate_ectopic_recombination(t, s)
    elseif event == IMMUNITY_LOSS
        get_rate_immunity_loss(t, s)
    elseif event == GENERALIZED_IMMUNITY_LOSS
        get_rate_generalized_immunity_loss(t, s)
    end
end

function do_event!(t, s, stats, event, db)
    if event == BITING
        do_biting!(t, s, stats)
    elseif event == IMMIGRATION
        do_immigration!(t, s, stats)
    elseif event == BACKGROUND_CLEARANCE
        do_background_clearance(t, s, stats)
    elseif event == SWITCHING
        do_switching!(t, s, stats)
    elseif event == MUTATION
        do_mutation!(t, s, stats)
    elseif event == ECTOPIC_RECOMBINATION
        do_ectopic_recombination!(t, s, stats)
    elseif event == IMMUNITY_LOSS
        do_immunity_loss!(t, s, stats)
    elseif event == GENERALIZED_IMMUNITY_LOSS
        do_generalized_immunity_loss!(t, s, stats)
    end
end


### BITING EVENT ###

function get_rate_biting(t, s)
    day_index = 1 + Int(floor(t)) % P.t_year
    biting_rate = if USE_BITING_RATE_MULTIPLIER_BY_YEAR
        year_index = 1 + Int(floor(t / P.t_year))
        P.biting_rate_multiplier_by_year[year_index] * P.biting_rate[day_index]
    else
        P.biting_rate[day_index]
    end
    biting_rate * P.n_hosts
end

function do_biting!(t, s, stats)

    stats.n_bites += 1
    s.n_bites_for_migration_rate += 1

    # Uniformly randomly sample infecting host (source) and host being infected
    # (destination).
    src_host = rand(s.hosts)
    dst_host = rand(s.hosts)

    # Advance host (rebirth or infection activation).
    advance_host!(t, s, src_host)
    advance_host!(t, s, dst_host)

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

    # Compute probability of each transmission.
    # Probability based on the generalized immunity level of the source host
    # and on a parameter which control how fast the GI impact the transmission.
    println("do_biting!: Compute probability of each transmission")
    transmitted_strains = []
    for active_infection in src_host.active_infections
        if P.generalized_immunity && src_host.generalized_immunity > 0
            if P.coinfection_reduces_transmission
                active_infection.p_transmit = P.transmissibility * exp(-1 * src_host.generalized_immunity * P.generalized_immunity_transmissibility) / src_active_count
            else
                active_infection.p_transmit = P.transmissibility * exp(-1 * src_host.generalized_immunity * P.generalized_immunity_transmissibility)
            end
            # Cost of carrying a resistance allele on the strain transmissibility.
            if P.drug_treatment
                if P.resistant_snp && active_infection.snps[1] == 2
                    active_infection.p_transmit = active_infection.p_transmit * P.resistant_cost
                end
            end
            println("The source host $(src_host.id) has $(src_active_count) active infections.")
            println("The source host $(src_host.id) has $(src_host.n_cleared_infections) cleared infections.")
            println("The source host $(src_host.id) has a GI level of $(src_host.generalized_immunity).")
            println("Its active infection $(active_infection.id) has a probability of transmission of $(active_infection.p_transmit).")
        else
            if P.coinfection_reduces_transmission
                active_infection.p_transmit = P.transmissibility / src_active_count
            else
                active_infection.p_transmit = P.transmissibility
            end
            # Cost of carrying a resistance allele on the strain transmissibility.
            if P.drug_treatment
                if P.resistant_snp && active_infection.snps[1] == 2
                    active_infection.p_transmit = active_infection.p_transmit * P.resistant_cost
                end
            end
            println("The source host $(src_host.id) has $(src_active_count) active infections.")
            println("The source host $(src_host.id) has $(src_host.n_cleared_infections) cleared infections.")
            println("Its active infection $(active_infection.id) has a probability of transmission of $(active_infection.p_transmit).")
        end

        # Choose if this active strain from the host will be transmitted to the mosquito.
        # This is determined by its transmissibility.
        if rand() < active_infection.p_transmit
            println("Its active infection $(active_infection.id) will be transmitted to the mosquito.")
            push!(transmitted_strains, active_infection)
        end
    end
    println("Total of $(length(transmitted_strains)) transmitted strains.")

    # The number of transmissions is bounded by the number of source infections
    # and the number of available slots in the destination.
    n_transmissions_max = min(length(transmitted_strains), dst_available_count)
    transmitted = false
    for i in 1:n_transmissions_max

        stats.n_transmissions += 1
        transmitted = true

        # Randomly sample two source infections within transmitted_strains to recombine.
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

        # Construct strain for new infection.
        if src_inf_1.strain_id == src_inf_2.strain_id
            # If both infections have the same strain, then the new infection
            # is given infection 1's genes with expression order shuffled.
            dst_inf.strain_id = src_inf_1.strain_id
            shuffle_columns_to!(dst_inf.genes, src_inf_1.genes)
            # The new infection is given infection 1's SNP alleles.
            if P.n_snps_per_strain > 0
                dst_inf.snps = src_inf_1.snps
            end
        else
            # Otherwise, the new infection is given a new strain constructed by
            # taking a random sample of the genes in the two source infections.
            dst_inf.strain_id = next_strain_id!(s)
            sample_columns_from_two_matrices_to!(dst_inf.genes, src_inf_1.genes, src_inf_2.genes)
            # The new infection is given a new set of SNP alleles constructed by
            # taking a random allele per SNP in the two source infections.
            if P.n_snps_per_strain > 0
                if !P.snp_linkage_disequilibrium
                    for i in 1:P.n_snps_per_strain
                        dst_inf.snps[i] = rand((src_inf_1.snps[i], src_inf_2.snps[i]))
                    end
                else
                    unlinked_snps = collect(1:P.n_snps_per_strain)
                    i = 1
                    while i in unlinked_snps

                        # Define the linked SNPs:
                        linked_snps = find_linked_snps(i)

                        # The new infection is given a new set of linked SNP alleles
                        # constructed by taking one allele per linked SNP in only
                        # one source infection.
                        if size(linked_snps)[1] > 1
                            unlinked_snps = setdiff(unlinked_snps, linked_snps)
                            parent_inf = rand((src_inf_1, src_inf_2))
                            for linked_snp in linked_snps
                                dst_inf.snps[linked_snp] = parent_inf.snps[linked_snp]
                            end
                        end
                        i += 1
                        while !(i in unlinked_snps) && i <= P.n_snps_per_strain
                            i += 1
                        end
                    end

                    # The new infection is given a new set of unlinked SNP alleles
                    # constructed by taking a random allele per unlinked SNP in
                    # the two source infections.
                    for unlinked_snp in unlinked_snps
                        dst_inf.snps[unlinked_snp] = rand((src_inf_1.snps[unlinked_snp], src_inf_2.snps[unlinked_snp]))
                    end
                end
            end
        end

        # For the new infection, calculate its detectability and pathogenicity.
        # Probabilities based on the generalized immunity level of its destination host
        # and on parameters controlling how fast the GI impact the detectability and pathogenicity.
        if P.generalized_immunity
            if dst_host.generalized_immunity > 0
                dst_inf.p_detect = P.detectability * exp(-1 * dst_host.generalized_immunity * P.generalized_immunity_detectability)
                println("The destination infection $(dst_inf.id) has a probability of detection of ($(dst_inf.p_detect)).")
                if P.drug_treatment
                    dst_inf.p_symptoms = P.pathogenicity * exp(-1 * dst_host.generalized_immunity * P.generalized_immunity_pathogenicity)
                    println("The destination infection $(dst_inf.id) has a probability of symptoms of ($(dst_inf.p_symptoms)).")
                    # Cost of carrying a resistance allele on detectability and pathogenicity.
                    if P.resistant_snp && dst_inf.snps[1] == 2
                        println("The destination infection $(dst_inf.id) carries the resistant allele ($(dst_inf.snps[1])).")
                        dst_inf.p_detect = dst_inf.p_detect * P.resistant_cost
                        dst_inf.p_symptoms = dst_inf.p_symptoms * P.resistant_cost
                        println("Its updated probability of detection is $(dst_inf.p_detect) and its probability of symptoms is $(dst_inf.p_symptoms).")
                    end
                end
            else
                # Note: make a function as same code used twice (see below).
                dst_inf.p_detect = P.detectability
                if P.drug_treatment
                    dst_inf.p_symptoms = P.pathogenicity
                    # Cost of carrying a resistance allele on detectability and pathogenicity.
                    if P.resistant_snp && dst_inf.snps[1] == 2
                        dst_inf.p_detect = dst_inf.p_detect * P.resistant_cost
                        dst_inf.p_symptoms = dst_inf.p_symptoms * P.resistant_cost
                    end
                end
            end
        else
            dst_inf.p_detect = P.detectability
            println("The destination infection $(dst_inf.id) has a probability of detection of ($(dst_inf.p_detect)).")
            if P.drug_treatment
                dst_inf.p_symptoms = P.pathogenicity
                println("The destination infection $(dst_inf.id) has a probability of symptoms of ($(dst_inf.p_symptoms)).")
                # Cost of carrying a resistance allele on detectability and pathogenicity.
                if P.resistant_snp && dst_inf.snps[1] == 2
                    println("The destination infection $(dst_inf.id) carries the resistant allele ($(dst_inf.snps[1])).")
                    dst_inf.p_detect = dst_inf.p_detect * P.resistant_cost
                    dst_inf.p_symptoms = dst_inf.p_symptoms * P.resistant_cost
                    println("Its updated probability of detection is $(dst_inf.p_detect) and its probability of symptoms is $(dst_inf.p_symptoms).")
                end
            end
        end
        # Add this infection to the destination host
        push!(dst_host.liver_infections, dst_inf)
        # Update population wide host liver max.
        s.n_liver_infections_per_host_max = max(s.n_liver_infections_per_host_max, length(dst_host.liver_infections))
    end

    if transmitted
        stats.n_transmitting_bites += 1
        s.n_transmitting_bites_for_migration_rate += 1
        true
    else
        false
    end
end

function advance_host!(t, s, host)
    if t > host.t_death
        # If the host is past its death time.
        do_rebirth!(t, s, host)
    else
        i = 1
        while i <= length(host.liver_infections)
            infection = host.liver_infections[i]
            if infection.t_infection + P.t_liver_stage < t

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

                    if P.generalized_immunity
                        # If detectability and transmissibility are affected by the generalized immunity (GI) level.
                        println("advance_host: GI is true")
                        if P.drug_treatment && rand() < infection.p_symptoms
                            # If infections are impacted by malaria drug treatments and this active infection is symptomatic.
                            # This choice is based on the probability that strain generates symptoms in its host, i.e. pathogenicity.
                            println("Infections are impacted by treatments and this infection $(infection.id), with a probability of generating symptoms of $(infection.p_symptoms), is symptomatic.")
                            if P.resistant_snp
                                # If there is one SNP under selection carrying a susceptible (1) or resistant (2) allele.
                                # Note: as it is used twice (i.e. for w/wo GI), we need to make a function!
                                if infection.snps[1] == 1
                                    # If this active infection carries a susceptible allele (1).
                                    # This active infections, and the other coinfecting strain(s) carrying a susceptible allele (1),
                                    # will be cleared without increasing the GI level. The coinfecting strain(s) carrying a resistant allele (2) will pursue.
                                    println("The infection $(infection.id), with allele $(infection.snps[1]) at the selected locus, is susceptible.")
                                    println("This infection is infecting host $(host.id).")
                                    println("Before drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                    println("Before drug treatment, host $(host.id) has a GI level of $(host.generalized_immunity).")
                                    if length(host.active_infections) > 1
                                        j = 1
                                        for active_infection in host.active_infections
                                            println("The coinfecting strain $(active_infection.id) has allele $(active_infection.snps[1]) at the selected locus.")
                                            if active_infection.snps[1] == 1
                                                println("The coinfecting strain $(active_infection.id) with allele $(active_infection.snps[1]) at the selected locus is susceptible.")
                                                delete_and_swap_with_end!(host.active_infections, j)
                                            end
                                            j += 1
                                        end
                                        k = 1
                                        for active_infection_detectable in host.active_infections_detectable
                                            if active_infection_detectable.snps[1] == 1
                                                delete_and_swap_with_end!(host.active_infections_detectable, k)
                                            end
                                            k += 1
                                        end
                                        println("After drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                        println("After drug treatment, host $(host.id) has a GI level of $(host.generalized_immunity).")
                                    else
                                        empty!(host.active_infections)
                                        empty!(host.active_infections_detectable)
                                        println("After drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                        println("After drug treatment, host $(host.id) has a GI level of $(host.generalized_immunity).")
                                    end
                                else
                                    # If this active infection carries a resistant allele (2).
                                    # Choose if this active strain could be detected.
                                    # This choice is based on its detectability in its host.
                                    println("The infection $(infection.id), with allele $(infection.snps[1]) at the selected locus, is resistant.")
                                    if rand() < infection.p_detect
                                        println("The infection $(infection.id), with a probability of detection of $(infection.p_detect), is detectable.")
                                        println("This infection is infecting host $(host.id) which currently has a GI level of $(host.generalized_immunity).")
                                        push!(host.active_infections_detectable, infection)
                                    else
                                        println("The infection $(infection.id), with a probability of detection of $(infection.p_detect), is undetectable.")
                                        println("This infection is infecting host $(host.id) which currently has a GI level of $(host.generalized_immunity).")
                                    end
                                    advance_immune_genes!(t, s, host, length(host.active_infections))
                                    if length(host.active_infections) > s.n_active_infections_per_host_max
                                        s.n_active_infections_per_host_max = length(host.active_infections)
                                        println("nb host.active_infections: $(length(host.active_infections))")
                                        println("s.n_active_infections_per_host_max: $(s.n_active_infections_per_host_max)")
                                    end
                                    # This active infections, and the other coinfecting strain(s) carrying a resistant allele (2) will pursue.
                                    # However, the coinfecting strain(s) carrying a susceptible allele (1) will be cleared without increasing the GI level.
                                    if length(host.active_infections) > 1
                                        j = 1
                                        for active_infection in host.active_infections
                                            println("The coinfecting strain $(active_infection.id) has allele $(active_infection.snps[1]) at the selected locus")
                                            if active_infection.snps[1] == 1
                                                println("The coinfecting strain $(active_infection.id) with allele $(active_infection.snps[1]) at the selected locus is susceptible")
                                                delete_and_swap_with_end!(host.active_infections, j)
                                            end
                                            j += 1
                                        end
                                        k = 1
                                        for active_infection_detectable in host.active_infections_detectable
                                            if active_infection_detectable.snps[1] == 1
                                                delete_and_swap_with_end!(host.active_infections_detectable, k)
                                            end
                                            k += 1
                                        end
                                    end
                                end
                            else
                                # If all the SNPs are neutral.
                                # All the active infections will be cleared without increasing the GI level.
                                println("All SNPs are neutral, i.e. all the active infections will be cleared.")
                                println("Before drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                println("Before drug treatment, host $(host.id) has a GI level of $(host.generalized_immunity).")
                                empty!(host.active_infections)
                                empty!(host.active_infections_detectable)
                                println("After drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                println("After drug treatment, host $(host.id) has a GI level of $(host.generalized_immunity).")
                            end
                        else
                            # If infections are not impacted by malaria drug treatments, or if this active infection is asymptomatic.
                            # Choose if the asymptomatic active strain could be detected.
                            # This choice is based on its detectability in its host.
                            println("Infections are not impacted by treatments, or this infection $(infection.id) is asymptomatic.")
                            println("This infection is infecting host $(host.id) which currently has a GI level of $(host.generalized_immunity).")
                            if rand() < infection.p_detect
                                println("This infection $(infection.id), with a probability of detection of $(infection.p_detect), is detectable.")
                                println("This infection is infecting host $(host.id) which currently has $(host.n_cleared_infections) cleared infections.")
                                push!(host.active_infections_detectable, infection)
                            else
                                println("This infection $(infection.id), with a probability of detection of $(infection.p_detect), is undetectable.")
                                println("This infection is infecting host $(host.id) which currently has $(host.n_cleared_infections) cleared infections.")
                            end
                            advance_immune_genes!(t, s, host, length(host.active_infections))
                            if length(host.active_infections) > s.n_active_infections_per_host_max
                                s.n_active_infections_per_host_max = length(host.active_infections)
                                println("nb host.active_infections: $(length(host.active_infections))")
                                println("s.n_active_infections_per_host_max: $(s.n_active_infections_per_host_max)")
                            end
                        end
                    else
                        # If detectability and transmissibility are not affected by the generalized immunity (GI) level.
                        println("advance_host: GI is false")
                        if P.drug_treatment && rand() < infection.p_symptoms
                            # If infections are impacted by malaria drug treatments and this active infection is symptomatic.
                            # This choice is based on the probability that strain generates symptoms in its host, i.e. pathogenicity.
                            println("Infections are impacted by treatments and this infection $(infection.id), with a probability of generating symptoms of $(infection.p_symptoms), is symptomatic.")
                            if P.resistant_snp
                                # If there is one SNP under selection carrying a susceptible (1) or resistant (2) allele.
                                # Note: as it is used twice (i.e. for w/wo GI), we need to make a function!
                                if infection.snps[1] == 1
                                    # If this active infection carries a susceptible allele (1).
                                    # This active infections, and the other coinfecting strain(s) carrying a susceptible allele (1),
                                    # will be cleared without increasing the GI level. The coinfecting strain(s) carrying a resistant allele (2) will pursue.
                                    println("The infection $(infection.id), with allele $(infection.snps[1]) at the selected locus, is susceptible.")
                                    println("This infection is infecting host $(host.id).")
                                    println("Before drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                    println("Before drug treatment, host $(host.id) has $(host.n_cleared_infections) cleared infections.")
                                    if length(host.active_infections) > 1
                                        j = 1
                                        for active_infection in host.active_infections
                                            println("The coinfecting strain $(active_infection.id) has allele $(active_infection.snps[1]) at the selected locus.")
                                            if active_infection.snps[1] == 1
                                                println("The coinfecting strain $(active_infection.id) with allele $(active_infection.snps[1]) at the selected locus is susceptible.")
                                                delete_and_swap_with_end!(host.active_infections, j)
                                            end
                                            j += 1
                                        end
                                        k = 1
                                        for active_infection_detectable in host.active_infections_detectable
                                            if active_infection_detectable.snps[1] == 1
                                                delete_and_swap_with_end!(host.active_infections_detectable, k)
                                            end
                                            k += 1
                                        end
                                        println("After drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                        println("After drug treatment, host $(host.id) has $(host.n_cleared_infections) cleared infections.")
                                    else
                                        empty!(host.active_infections)
                                        empty!(host.active_infections_detectable)
                                        println("After drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                        println("After drug treatment, host $(host.id) has $(host.n_cleared_infections) cleared infections.")
                                    end
                                else
                                    # If this active infection carries a resistant allele (2).
                                    # Choose if this active strain could be detected.
                                    # This choice is based on its detectability in its host.
                                    println("The infection $(infection.id), with allele $(infection.snps[1]) at the selected locus, is resistant.")
                                    if rand() < infection.p_detect
                                        println("The infection $(infection.id), with a probability of detection of $(infection.p_detect), is detectable.")
                                        println("This infection is infecting host $(host.id) which currently has $(host.n_cleared_infections) cleared infections.")
                                        push!(host.active_infections_detectable, infection)
                                    else
                                        println("The infection $(infection.id), with a probability of detection of $(infection.p_detect), is undetectable.")
                                        println("This infection is infecting host $(host.id) which currently has $(host.n_cleared_infections) cleared infections.")
                                    end
                                    advance_immune_genes!(t, s, host, length(host.active_infections))
                                    if length(host.active_infections) > s.n_active_infections_per_host_max
                                        s.n_active_infections_per_host_max = length(host.active_infections)
                                        println("nb host.active_infections: $(length(host.active_infections))")
                                        println("s.n_active_infections_per_host_max: $(s.n_active_infections_per_host_max)")
                                    end
                                    # This active infections, and the other coinfecting strain(s) carrying a resistant allele (2) will pursue.
                                    # However, the coinfecting strain(s) carrying a susceptible allele (1) will be cleared without increasing the GI level.
                                    if length(host.active_infections) > 1
                                        j = 1
                                        for active_infection in host.active_infections
                                            println("The coinfecting strain $(active_infection.id) has allele $(active_infection.snps[1]) at the selected locus.")
                                            if active_infection.snps[1] == 1
                                                println("The coinfecting strain $(active_infection.id), with allele $(active_infection.snps[1]) at the selected locus, is susceptible.")
                                                delete_and_swap_with_end!(host.active_infections, j)
                                            end
                                            j += 1
                                        end
                                        k = 1
                                        for active_infection_detectable in host.active_infections_detectable
                                            if active_infection_detectable.snps[1] == 1
                                                delete_and_swap_with_end!(host.active_infections_detectable, k)
                                            end
                                            k += 1
                                        end
                                    end
                                end
                            else
                                # If all the SNPs are neutral.
                                # All the active infections will be cleared without increasing the GI level.
                                println("All SNPs are neutral, i.e. all the active infections will be cleared.")
                                println("Before drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                #println("Before drug treatment, host $(host.id) has a GI level of $(host.generalized_immunity).")
                                empty!(host.active_infections)
                                empty!(host.active_infections_detectable)
                                println("After drug treatment, host $(host.id) has $(length(host.active_infections)) active infections.")
                                #println("After drug treatment, host $(host.id) has a GI level of $(host.generalized_immunity).")
                            end
                        else
                            # If infections are not impacted by malaria drug treatments, or if this active infection is asymptomatic.
                            # Choose if the asymptomatic active strain could be detected.
                            # This choice is based on its detectability in its host.
                            println("Infections are not impacted by treatments, or the infection $(infection.id) is asymptomatic.")
                            if rand() < infection.p_detect
                                println("This infection $(infection.id), with a probability of detection of $(infection.p_detect), is detectable.")
                                println("This infection is infecting host $(host.id) which currently has $(host.n_cleared_infections) cleared infections.")
                                push!(host.active_infections_detectable, infection)
                            else
                                println("This infection $(infection.id), with a probability of detection of $(infection.p_detect), is undetectable.")
                                println("This infection is infecting host $(host.id) which currently has $(host.n_cleared_infections) cleared infections.")
                            end
                            advance_immune_genes!(t, s, host, length(host.active_infections))
                            if length(host.active_infections) > s.n_active_infections_per_host_max
                                s.n_active_infections_per_host_max = length(host.active_infections)
                                println("nb host.active_infections: $(length(host.active_infections))")
                                println("s.n_active_infections_per_host_max: $(s.n_active_infections_per_host_max)")
                            end
                        end
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
    empty!(host.active_infections_detectable)
    empty!(host.immunity)
    host.generalized_immunity = 0
end


### IMMIGRATION EVENT ###

function get_rate_immigration(t, s)
    P.immigration_rate_fraction * get_rate_biting(t, s) * s.infected_ratio
end

function do_immigration!(t, s, stats)

    # Sample a random host and advance it (rebirth or infection activation).
    host = rand(s.hosts)
    advance_host!(t, s, host)

    # If host doesn't have an available infection slot, reject this sample.
    if !isnothing(P.n_infections_liver_max)
        if length(host.liver_infections) == P.n_infections_liver_max
            return false
        end
    end

    # Construct infection by sampling from gene pool.
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
    for snp in 1:P.n_snps_per_strain
        infection.snps[snp] = sample([1, 2],
        Weights([s.initial_snp_allele_frequencies[snp], 1 - s.initial_snp_allele_frequencies[snp]]))
    end

    # Add infection to host.
    push!(host.liver_infections, infection)

    s.n_liver_infections_per_host_max = max(s.n_liver_infections_per_host_max, length(host.liver_infections))

    true
end


### RANDOM BACKGROUND CLEARANCE EVENT ###

function get_rate_background_clearance(t, s)
    # The total rate includes both active and liver infections because host state may not be fully up to date,
    # and a liver infection may be activated when host state is updated to the current time.
    # Rejection sampling is used to effect the correct rate.
    P.background_clearance_rate * P.n_hosts * (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max)
end

function do_background_clearance(t, s, stats)
    index = rand(CartesianIndices((P.n_hosts, (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max))))
    host = s.hosts[index[1]]

    # Advance host (rebirth or infection activation).
    advance_host!(t, s, host)
    inf_index = index[2]

    # If the infection index is out of range, this is a rejected sample.
    # Otherwise we'll proceed.
    if inf_index > length(host.active_infections)
        false
    else
        infection = host.active_infections[inf_index]
        clear_active_infection!(t, s, host, inf_index)
        inf_det_index = 1
        for inf_det in host.active_infections_detectable
            println("One of the detectable infection: $(inf_det.id)")
            if infection.id == inf_det.id
                println("Infection $(inf_det.id) corresponds to the active infection!")
                act_inf_det = host.active_infections_detectable[inf_det_index]
                println("Confirmation that infection $(act_inf_det.id) index is $(inf_det_index).")
                delete_and_swap_with_end!(host.active_infections_detectable, inf_det_index)
            end
            inf_det_index += 1
        end

        true
    end
end


### SWITCHING EVENT ###

function get_rate_switching(t, s)
    # The total rate includes both active and liver infections because host state may not be fully up to date,
    # and a liver infection may be activated when host state is updated to the current time.
    # Rejection sampling is used to effect the correct rate.
    if !P.whole_gene_immune
        # Switching rate set by total number of alleles.
        (P.switching_rate * P.n_loci) * P.n_hosts * (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max)
    else
        P.switching_rate * P.n_hosts * (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max)
    end
end

function do_switching!(t, s, stats)
    println("do_switching! actually happening, n_hosts: $(P.n_hosts)")
    println("do_switching! actually happening, n_active_infections_per_host_max: $(s.n_active_infections_per_host_max)")
    println("do_switching! actually happening, n_liver_infections_per_host_max: $(s.n_liver_infections_per_host_max)")
    index = rand(CartesianIndices((P.n_hosts, (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max))))
    host = s.hosts[index[1]]

    # Advance host (rebirth or infection activation).
    advance_host!(t, s, host)
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

    # If we're at the end, clear the infection and return.
    if infection.expression_index == P.n_genes_per_strain && infection.expression_index_locus == P.n_loci
        clear_active_infection!(t, s, host, inf_index)
        # Make another "clear_active_infection" function (but with GI).
        host.generalized_immunity += 1
        inf_det_index = 1
        for inf_det in host.active_infections_detectable
            println("One of the detectable infection: $(inf_det.id)")
            if infection.id == inf_det.id
                println("Infection $(inf_det.id) corresponds to the active infection!")
                act_inf_det = host.active_infections_detectable[inf_det_index]
                println("Confirmation that infection $(act_inf_det.id) index is $(inf_det_index).")
                delete_and_swap_with_end!(host.active_infections_detectable, inf_det_index)
            end
            inf_det_index += 1
        end
    else
        # Otherwise, advance gene and/or locus expression(s).
        if !P.whole_gene_immune
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

    """
    after gaining immunity to the expressed gene, loop through the other infections
    in the host to see if any one needs advancing
    """
    i = 1
    while i <= length(host.active_infections)
        if advance_immune_genes!(t,s,host,i)
            # If there is no end of expression and reordering of infections, then index plus 1.
            i += 1
        end
    end

    return true
end

# This function moves the expression index of an infection to its first non-immune allele/gene.
function advance_immune_genes!(t, s, host, inf_index)
    # Advance expression until a non-immune gene or allele is reached.
    infection = host.active_infections[inf_index]
    # If the host not immune, stop advancing.
    to_advance = is_immune(host.immunity, infection.genes[:, infection.expression_index],infection.expression_index_locus)

    while to_advance
        # Increment immunity level to currently expressed gene or allele.
        if infection.expression_index_locus == P.n_loci
            increment_immunity!(s, host, infection.genes[:, infection.expression_index])
        end

        # If we're at the end, clear the infection and return.
        if infection.expression_index == P.n_genes_per_strain && infection.expression_index_locus == P.n_loci
            clear_active_infection!(t, s, host, inf_index)
            host.generalized_immunity += 1
            inf_det_index = 1
            for inf_det in host.active_infections_detectable
                println("One of the detectable infection: $(inf_det.id)")
                if infection.id == inf_det.id
                    println("Infection $(inf_det.id) corresponds to the active infection!")
                    act_inf_det = host.active_infections_detectable[inf_det_index]
                    println("Confirmation that infection $(act_inf_det.id) index is $(inf_det_index).")
                    delete_and_swap_with_end!(host.active_infections_detectable, inf_det_index)
                end
                inf_det_index += 1
            end
            # Here if deleting an infection, and the indexing changed, then tell the calling function
            # it doesn't advancing its index.
            return false
        else
            # Otherwise, advance gene and/or locus expression(s).
            if !P.whole_gene_immune
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
# Update mutation and recombination rates towards all infections.
function get_rate_mutation(t, s)
    # The total rate includes both active and liver infections because host state may not be fully up to date,
    # and a liver infection may be activated when host state is updated to the current time.
    # Rejection sampling is used to effect the correct rate.
    P.mutation_rate * P.n_hosts * (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max) * P.n_genes_per_strain * P.n_loci
end

function do_mutation!(t, s, stats)
    index = rand(CartesianIndices((P.n_hosts, (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max), P.n_genes_per_strain, P.n_loci)))
    host = s.hosts[index[1]]
    inf_index = index[2]
    expression_index = index[3]
    locus = index[4]

    # Advance host (rebirth or infection activation).
    advance_host!(t, s, host)

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
    s.next_gene_id_mut += 1

    true
end


### ECTOPIC RECOMBINATION EVENT ###

function get_rate_ectopic_recombination(t, s)
    # The total rate includes both active and liver infections because host state may not be fully up to date,
    # and a liver infection may be activated when host state is updated to the current time.
    # Rejection sampling is used to effect the correct rate.
    P.ectopic_recombination_rate *
        P.n_hosts * (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max) *
        P.n_genes_per_strain * (P.n_genes_per_strain - 1) / 2.0
end

function do_ectopic_recombination!(t, s, stats)
    # Index based on the total number of infections.
    index = rand(CartesianIndices((P.n_hosts, (s.n_active_infections_per_host_max + s.n_liver_infections_per_host_max))))
    host = s.hosts[index[1]]
    inf_index = index[2]

    # Advance host (rebirth or infection activation).
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

    # If the genes are the same, this is a no-op.
    if gene1 == gene2
        return false
    end

    breakpoint, p_functional = if P.ectopic_recombination_generates_new_alleles
        # Choose a breakpoint.
        breakpoint = P.n_loci * rand()
        p_functional = p_recombination_is_functional_real(gene1, gene2, breakpoint)

        (Int(ceil(breakpoint)), p_functional)
    else
        # Choose a breakpoint.
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

    # Recombine to modify first gene, if functional.
    if !is_conversion && rand() < p_functional
        infection.genes[:, gene_index_1] = if create_new_allele
            recombine_genes_new_allele(s, gene1, gene2, breakpoint)
        else
            recombine_genes(gene1, gene2, breakpoint, s)
        end
    end

    # Recombine to modify second gene, if functional.
    if rand() < p_functional
        infection.genes[:, gene_index_2] = if create_new_allele
            recombine_genes_new_allele(s, gene2, gene1, breakpoint)
        else
            recombine_genes(gene2, gene1, breakpoint, s)
        end
    end

    if gene1 == infection.genes[:, gene_index_1] && gene2 == infection.genes[:, gene_index_2]
        recombined = false
    else
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
    Clear active infection
    This function is used by all code where active infections can be cleared:
    * Clearance at end of expression sequence during a switching event
    * Clearance
    * Random background clearance events
"""
function clear_active_infection!(t, s, host, inf_index)
    if s.n_cleared_infections % P.sample_infection_duration_every == 0
        # Calculate and write the infection duration.
        add_infection_duration!(t, s, host, inf_index)
    end
    s.n_cleared_infections += 1
    host.n_cleared_infections += 1
    push!(s.old_infections, host.active_infections[inf_index])
    delete_and_swap_with_end!(host.active_infections, inf_index)
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

function recombine_genes(gene1, gene2, breakpoint, s)
    gene = MGene(undef)
    gene[1:(breakpoint - 1)] = gene1[1:(breakpoint - 1)]
    gene[breakpoint:end] = gene2[breakpoint:end]
    if gene != gene1 && gene != gene2
        s.next_gene_id_recomb += 1
    end
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
    if gene != gene1 && gene != gene2
        s.next_gene_id_recomb += 1
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

    # Advance host (rebirth or infection activation).
    advance_host!(t, s, host)

    # If the immunity index is beyond this host's immunity count, reject this sample.
    if immunity_index >  immunity_count(host.immunity)
        return false
    end

    decrement_immunity_at_sampled_index!(host.immunity, immunity_index)

    false
end


### GENERALIZED IMMUNITY LOSS EVENT ###

function get_rate_generalized_immunity_loss(t, s)
    P.generalized_immunity_loss_rate * P.n_hosts
end

function do_generalized_immunity_loss!(t, s, stats)
    println("Do generalized immunity loss.")
    index = rand(1:P.n_hosts)
    host = s.hosts[index]
    println("Host involved in the generalized immunity loss: $(host.id)")

    # Advance host (rebirth or infection activation).
    advance_host!(t, s, host)

    # If the host generalized immunity level is null, reject this sample.
    if host.generalized_immunity < 1
        println("The host $(host.id) GI level is $(host.generalized_immunity) and cannot be decremented.")
        return false
    end

    # Decrement generalized immunity.
    println("The host $(host.id) GI level is $(host.generalized_immunity) and can be decremented.")
    println("The host $(host.id) has $(host.n_cleared_infections) cleared infections.")
    host.generalized_immunity -= 1
    println("The host $(host.id) GI level after loss is $(host.generalized_immunity).")
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

"""
    Find the SNPs that are in linkage disequilibrium (LD) with SNP i (i.e. linked SNPs).
    For each SNP, the function selects the linked SNP using the LD coefficients
    from the pairwise LD matrix to weight the probability to drawing a given SNP.
"""
function find_linked_snps(i)
    linked_snps = [i]
    for j = (i+1):P.n_snps_per_strain
        if rand() < P.snp_pairwise_ld[j, i]
            append!(linked_snps, j)
        end
    end
    linked_snps
end


### IMMUNITY FUNCTIONS ###

function empty!(ih::ImmuneHistory)
    for d in ih.vd
        empty!(d)
    end
end

function immunity_count(ih::ImmuneHistory)
    sum(length(ih.vd[locus]) for locus in 1:length(ih.vd))
end


function increment_immunity!(s, host, gene)
    increment_immunity!(host.immunity, gene)
    s.n_immunities_per_host_max = max(s.n_immunities_per_host_max, immunity_count(host.immunity))
end

function increment_immunity!(ih::ImmuneHistory, gene)
    # Increment immunity at each locus.
    for (locus, allele_id) in enumerate(gene)
        old_level = get(ih.vd[locus], allele_id, ImmunityLevel(0))
        if old_level < typemax(ImmunityLevel)
            ih.vd[locus][allele_id] = old_level + 1
        end
    end
end

function decrement_immunity_at_sampled_index!(ih::ImmuneHistory, index)
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

function decrement_immunity!(ih::ImmuneHistory, locus::Locus, allele_id::AlleleId)
    # Get old immunity level from immunity dict.
    old_level = ih.vd[locus][allele_id]

    # Decrement immunity, removing it entirely if we reach 0.
    if old_level == 1
        delete!(ih.vd[locus], allele_id)
    else
        ih.vd[locus][allele_id] -= 1
    end
end

function is_immune(ih::ImmuneHistory, gene, loc)
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

# Push the current infection's duration calculations into durations vector.
function add_infection_duration!(t, s, host, i)
    get_duration!(host.active_infections, i, t)
    newInfDur = InfectionDuration(
        id=host.active_infections[i].id,
        host_id=host.id,
        n_cleared_infections=host.n_cleared_infections,
        n_immune_alleles=immunity_count(host.immunity),
        t_infection = host.active_infections[i].t_infection,
        t_expression = host.active_infections[i].t_expression,
        duration = host.active_infections[i].duration
    )
    push!(s.durations, newInfDur)
end
