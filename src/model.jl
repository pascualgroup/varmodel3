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

include("util.jl")
include("state.jl")
include("output.jl")

import Profile
import Serialization

const N_EVENTS = 9
const EVENTS = collect(1:N_EVENTS)
const (DEATH, BITING, IMMIGRATION, BACKGROUND_CLEARANCE, LIVER_PROGRESS, SWITCHING, MUTATION, ECTOPIC_RECOMBINATION, IMMUNITY_LOSS) = EVENTS

const USE_BITING_RATE_MULTIPLIER_BY_YEAR = P.biting_rate_multiplier_by_year !== nothing

gene_pair_indices = []
for i in 1:(P.n_genes_per_strain-1)
    for j in (i+1):P.n_genes_per_strain
        push!(gene_pair_indices, (i,j))
    end
end

num_genes_var_groups = []
for i in 1:length(P.var_groups_ratio)
    num_genes_var_group = round(Int, P.var_groups_ratio[i] * P.n_genes_initial) 
    push!(num_genes_var_groups, num_genes_var_group)
end

allele_ids_var_groups = []
num_allele_ids_var_groups = round.(Int, P.var_groups_ratio * P.n_alleles_per_locus_initial)
for i in 1:length(P.var_groups_ratio)
    allele_ids_var_group = if i == 1
        1:num_allele_ids_var_groups[i]
    else 
        (sum(num_allele_ids_var_groups[1:(i-1)])+1):sum(num_allele_ids_var_groups[1:i])
    end
    push!(allele_ids_var_groups, allele_ids_var_group)
end

infection_genes_index_var_groups = []
for i in 1:length(P.var_groups_ratio)
    infection_genes_index_var_group = if i == 1
        1:round(Int, P.var_groups_ratio[i] * P.n_genes_per_strain)
    else
        (sum(round.(Int, P.var_groups_ratio[1:(i-1)] * P.n_genes_per_strain)) + 1):sum(round.(Int, P.var_groups_ratio[1:i] * P.n_genes_per_strain))
    end
    push!(infection_genes_index_var_groups, infection_genes_index_var_group)
end

func_rank = ordinalrank(P.var_groups_functionality, rev = true)

function run()
    if P.profile_on
        profile()
    else
        run_inner()
    end
end

function profile()
    println("Running with profiling on...")
    Profile.init(n = 10^7, delay = P.profile_delay)
    @Profile.profile run_inner()

    println("About to write profile...")
    profile_data = Profile.retrieve()
    Serialization.serialize(P.profile_filename, profile_data)
    println("Profile written.")
end

function run_inner()
    db = initialize_database()

    # Seed the random number generator using the provided seed,
    # or, if absent, by generating one from the OS's source of entropy.
    rng_seed = if isnothing(P.rng_seed)
        rand(RandomDevice(), 1:typemax(Int64))
    else
        P.rng_seed
    end
    # Random.seed!(rng_seed)
    rng = Xoshiro(rng_seed)
    execute(db.meta, ("rng_seed", rng_seed))

    # Start recording elapsed time.
    start_datetime = now()

    # Used to decide when to do output and verification.
    t_next_integer_float = 1.0

    # Initialize state.
    t = 0.0
    s = initialize_state(rng)
    verify(t, s)

    # Run initial output.
    stats = SummaryStats()
    write_output!(db, 0, s, stats)

    # Initialize event rates.
    # total_rate = sum(rates)
    #weights = Weights(rates, total_rate)
    event_dist = WeightedDiscreteDistribution(10.0, [get_rate(t, s, event) for event in EVENTS])

    # Batched exponential distribution for event loop draws
    batched_exp_dist = BatchedDistribution(Exponential(1.0), P.rng_batch_size)

    # Loop events until end of simulation.
    while total_weight(event_dist) > 0.0 && t < P.t_end
        # Draw next time with rate equal to the sum of all event rates.
        dt = rand(rng, batched_exp_dist) / total_weight(event_dist)
        # @assert dt > 0.0 && !isinf(dt)

        # At each integer time, write output/state verification (if necessary),
        # and update the biting/immigration rate.
        # Loop required in case the simulation jumps past two integer times.
        t_next = t + dt
        while t_next_integer_float < t_next
            t_next_integer = Int64(t_next_integer_float)
            write_output!(db, t_next_integer, s, stats)
            if P.verification_period !== nothing && t_next_integer % P.verification_period == 0
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
            
            if t_next_integer % P.gene_group_id_association_recomputation_period == 0
                recompute_gene_group_id_association!(s)
            end

            # Update all rates & reset rate total to prevent error accumulation
            for event in EVENTS
                update_rate!(t_next_integer, s, event_dist, event)
            end
            recompute_total_weight!(event_dist)
            
            t_next_integer_float += 1.0
        end

        # Draw the event, update time, and execute event.
        # event = direct_sample_linear_scan(rates, total_rate)
        # event = sample(weights)
        event = rand(rng, event_dist)
        t = t_next
        if do_event!(t, s, stats, event, event_dist)
            stats.n_events += 1
        end
    end

    elapsed_time = Dates.value(now() - start_datetime) / 1000.0
    println("elapsed time (s): $(elapsed_time)")
    execute(db.meta, ("elapsed_time", elapsed_time))
    
    max_rss_gb = Sys.maxrss() / 2^30
    println("maxrss (GB) = $(max_rss_gb)")
    execute(db.meta, ("max_rss_gb", max_rss_gb))

    went_extinct = total_weight(event_dist) == 0.0
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


# Remove extinct genes from the association of gene to group_id dictionary.
function recompute_gene_group_id_association!(s)
    liver_infection_genes = Set()
    active_infection_genes = Set()
    for host in s.hosts
        for liver_infection in host.liver_infections
            for i in 1:size(liver_infection.genes)[2]
                gene_temp_alleles = liver_infection.genes[:,i]
                gene_temp = Gene(gene_temp_alleles)
                push!(liver_infection_genes, gene_temp)
            end
        end
        
        for active_infection in host.active_infections
            for i in 1:size(active_infection.genes)[2]
                gene_temp_alleles = active_infection.genes[:,i]
                gene_temp = Gene(gene_temp_alleles)
                push!(active_infection_genes, gene_temp)
            end
        end
    end
    pool_genes = Set()
    for i in 1:size(s.gene_pool)[2]
        gene_temp_alleles = s.gene_pool[:, i]
        gene_temp = Gene(gene_temp_alleles)
        push!(pool_genes, gene_temp)
    end
    genes_circulating = union(liver_infection_genes, active_infection_genes, pool_genes)
    genes_in_dictionary = collect(keys(s.association_genes_to_var_groups))
    genes_extinct = setdiff(genes_in_dictionary, genes_circulating)
    
    for gene_extinct in genes_extinct 
        delete!(s.association_genes_to_var_groups, gene_extinct)
    end
end


### INITIALIZATION ###

#function initialize_state(a)
function initialize_state(rng)
    println("initialize_state()")

    # Initialize gene pool as an (n_loci, n_genes_initial) matrix whose columns
    # are unique randomly generated genes. Initialize gene-to-group-id map as an empty dictionary.
    gene_pool_set = Set()
    association_genes_to_var_groups_init = Dict{Gene, GeneGroupId}()
    for group_id in 1:length(P.var_groups_ratio) # Create genes for each group
        num_genes_var_group = num_genes_var_groups[group_id] 
        if num_genes_var_group > 0
            counter = sum(num_genes_var_groups[1:group_id])
            if !P.var_groups_do_not_share_alleles 
                allele_ids_var_group = 1:P.n_alleles_per_locus_initial
            else  # If genes from different group do not share allele, then split alleles into groups corresponding to gene groups.
                allele_ids_var_group = allele_ids_var_groups[group_id]
            end
            
            while length(gene_pool_set) < counter
                gene_temp_alleles = rand(rng, allele_ids_var_group, P.n_loci)
                if !(gene_temp_alleles in gene_pool_set)
                    push!(gene_pool_set, gene_temp_alleles)
                    gene_temp = Gene(gene_temp_alleles)
                    association_genes_to_var_groups_init[gene_temp] = group_id
                    # @assert haskey(association_genes_to_var_groups_init, gene_temp)
                end
            end
        end
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
            t_birth = 0, # t_death = draw_host_lifetime(),
            liver_infections = [], active_infections = [],
            immunity = ImmuneHistory(),
            n_cleared_infections = 0
        )
        for id in 1:P.n_hosts
    ]
    for host in hosts
        # Chosen to roughly match the behavior of the old code, so some hosts die off "early"
        # in case we want to compare lifetime distributions,
        # but it doesn't actually affect dynamics at all since it's a Poisson process
        host.t_birth = -rand(rng) * rand(rng, Exponential(P.mean_host_lifetime))
    end

    # Infect n_initial_infections hosts at t = 0. Genes in the infection
    # are sampled uniformly randomly from the gene pool.
    # Check for var_groups_fix_ratio
    for (infection_id, host_index) in enumerate(sample(rng, 1:P.n_hosts, P.n_initial_infections, replace = false))
        infection = create_empty_infection()
        infection.id = infection_id
        infection.t_infection = 0.0
        infection.liver_index = 1
        infection.expression_index = 0
        infection.t_expression = NaN
        infection.duration = NaN
        infection.strain_id = infection_id
        if !P.var_groups_fix_ratio
            infection.genes[:,:] = reshape(
                gene_pool[:, rand(rng, 1:P.n_genes_initial, P.n_genes_per_strain)],
                (P.n_loci, P.n_genes_per_strain)
            )
            group_ids = []
            for i in 1:size(infection.genes)[2]
                group_id = association_genes_to_var_groups_init[Gene(infection.genes[:,i])]
                push!(group_ids, group_id)
            end
        else
            infection_genes = zeros(AlleleId, P.n_loci, P.n_genes_per_strain)
            group_ids = []
            for group_id in 1:length(P.var_groups_ratio)
                if P.var_groups_ratio[group_id] > 0.0
                    infection_genes_index_var_group = infection_genes_index_var_groups[group_id]
                    for infection_gene_index_var_group in infection_genes_index_var_group
                        infection_gene = gene_pool[:,rand(rng, 1:P.n_genes_initial, 1)]
                        infection_gene_group_id = association_genes_to_var_groups_init[Gene(infection_gene)]
                        while infection_gene_group_id != group_id # Keeping drawing until the targeted group of genes is drawn.
                            infection_gene = gene_pool[:,rand(rng, 1:P.n_genes_initial, 1)]
                            infection_gene_group_id = association_genes_to_var_groups_init[Gene(infection_gene)]
                        end
                        # @assert association_genes_to_var_groups_init[Gene(infection_gene)] == group_id
                        push!(group_ids, group_id)
                        infection_genes[:,infection_gene_index_var_group] = infection_gene # sample with replacement
                    end
                end
            end
            infection.genes[:,:] = reshape(infection_genes, (P.n_loci, P.n_genes_per_strain))
        end

        if P.var_groups_high_functionality_express_earlier
            genes_reorder = zeros(AlleleId, P.n_loci, P.n_genes_per_strain)
            index_temp_all = []
            for i in 1:length(func_rank)
                group_id_temp = findfirst(item->item == i, func_rank)
                index_temp = findall(item->item == group_id_temp, group_ids)
                push!(index_temp_all, length(index_temp))
                if length(index_temp) > 0
                    index_temp_reorder = if i == 1
                        1:length(index_temp)
                    else 
                        (sum(index_temp_all[1:(i-1)])+1):sum(index_temp_all[1:i])
                    end
                    genes_reorder[:,index_temp_reorder] = infection.genes[:, index_temp]
                end
            end
            infection.genes[:,:] = reshape(genes_reorder, (P.n_loci, P.n_genes_per_strain))
        end

        push!(hosts[host_index].liver_infections, infection)
    end

    State(
        rng = rng,
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
        association_genes_to_var_groups = association_genes_to_var_groups_init
    )
end

"""
    Draws host lifetime from a distribution.

    The distribution is an exponential distribution with mean
    `mean_host_lifetime`, truncated at `max_host_lifetime`.
"""
function draw_host_lifetime(rng)
    dist = Exponential(P.mean_host_lifetime)
    while true
        lifetime = rand(rng, dist)
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
        strain_id = StrainId(0),
        genes = fill(AlleleId(0), (P.n_loci, P.n_genes_per_strain)),
        liver_index = LiverIndex(0),
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

function update_rate!(t, s, event_dist, event)
    update!(event_dist, event, get_rate(t, s, event))
end

function update_rates_from_n_active_infections_per_host_max!(t, s, event_dist)
    update_rate!(t, s, event_dist, SWITCHING)
    update_rate!(t, s, event_dist, BACKGROUND_CLEARANCE)
    update_rate!(t, s, event_dist, MUTATION)
    update_rate!(t, s, event_dist, ECTOPIC_RECOMBINATION)
    update_rate!(t, s, event_dist, BACKGROUND_CLEARANCE)
end

function get_rate(t, s, event)
    if event == DEATH
        get_rate_death(t, s)
    elseif event == BITING
        get_rate_biting(t, s)
    elseif event == IMMIGRATION
        get_rate_immigration(t, s)
    elseif event == BACKGROUND_CLEARANCE
        get_rate_background_clearance(t, s)
    elseif event == LIVER_PROGRESS
        get_rate_liver_progress(t, s)
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

function do_event!(t, s, stats, event, event_dist)
    if event == DEATH
        do_death!(t, s, stats, event_dist)
    elseif event == BITING
        do_biting!(t, s, stats, event_dist)
    elseif event == IMMIGRATION
        do_immigration!(t, s, stats, event_dist)
    elseif event == BACKGROUND_CLEARANCE
        do_background_clearance(t, s, stats, event_dist)
    elseif event == LIVER_PROGRESS
        do_liver_progress!(t, s, stats, event_dist)
    elseif event == SWITCHING
        do_switching!(t, s, stats, event_dist)
    elseif event == MUTATION
        do_mutation!(t, s, stats, event_dist)
    elseif event == ECTOPIC_RECOMBINATION
        do_ectopic_recombination!(t, s, stats, event_dist)
    elseif event == IMMUNITY_LOSS
        do_immunity_loss!(t, s, stats, event_dist)
    end
end


### DEATH EVENT ###

function get_rate_death(t, s)
    P.n_hosts / P.mean_host_lifetime
end

function do_death!(t, s, stats, event_dist)
    host = rand(s.rng, s.hosts)
    do_rebirth!(t, s, host)
    true
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

function do_biting!(t, s, stats, event_dist)
    stats.n_bites += 1
    s.n_bites_for_migration_rate += 1

    # Uniformly randomly sample infecting host (source) and host being infected
    # (destination).
    src_host = rand(s.rng, s.hosts)
    dst_host = rand(s.rng, s.hosts)

    # The source host must be infected in order to transmit.
    src_active_count = length(src_host.active_infections)
    if src_active_count == 0
        return false
    end
    stats.n_infected_bites += 1
    s.n_transmitting_bites_for_migration_rate += 1

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
    p_transmit = if P.coinfection_reduces_transmission
        P.transmissibility / src_active_count
    else
        P.transmissibility
    end

    # First choose the active strains from the host that will be transmitted to mosquito.
    # This is determined by the transmissibility.
    choose_transmit = rand(s.rng, Float64, src_active_count)
    
    # Find out the currently expressed gene, and its group id, which impacts the transmissibility of the infection.
    infs_transmissibility = []
    for inf_temp in src_host.active_infections
        gene_index = inf_temp.expression_index
        # @assert 1 <= gene_index <= P.n_genes_per_strain
        gene_temp = Gene(inf_temp.genes[:,gene_index])
        # @assert haskey(s.association_genes_to_var_groups, gene_temp)
        gene_temp_group_id = s.association_genes_to_var_groups[gene_temp]
        inf_transmissibility = P.var_groups_functionality[gene_temp_group_id]
        push!(infs_transmissibility, inf_transmissibility)
    end
    
    transmitted_strains = src_host.active_infections[choose_transmit.<p_transmit*infs_transmissibility]
    #println("t = $(t): originalSize $(src_active_count), newSize $(length(transmitted_strains))")

    # The number of transmissions is bounded by the number of source infections
    # and the number of available slots in the destination.
    n_transmissions_max = min(length(transmitted_strains), dst_available_count)
    transmitted = false
    should_update_rates = false
    for i in 1:n_transmissions_max
        stats.n_transmissions += 1
        transmitted = true

        # Randomly sample two source infections within transmitted_strains to recombine.
        src_inf_1 = rand(s.rng, transmitted_strains)
        src_inf_2 = rand(s.rng, transmitted_strains)

        # Get a new infection struct, or recycle an old infection
        # to prevent excess memory allocation.
        dst_inf = recycle_or_create_infection(s)
        dst_inf.id = next_infection_id!(s)
        dst_inf.t_infection = t
        dst_inf.t_expression = NaN
        dst_inf.liver_index = 1
        dst_inf.expression_index = 0
        dst_inf.expression_index_locus = 0
        dst_inf.duration = NaN

        # Construct strain for new infection.
        if src_inf_1.strain_id == src_inf_2.strain_id
            # If both infections have the same strain, then the new infection
            # is given infection 1's genes with expression order shuffled.
            dst_inf.strain_id = src_inf_1.strain_id
            shuffle_columns_to!(s.rng, dst_inf.genes, src_inf_1.genes)
        else
            # Otherwise, the new infection is given a new strain constructed by
            # taking a random sample of the genes in the two source infections.
            dst_inf.strain_id = next_strain_id!(s)
            sample_columns_from_two_matrices_to_util2!(s.rng, dst_inf.genes, src_inf_1.genes, src_inf_2.genes, P, s, infection_genes_index_var_groups)
        end

        # Add this infection to the destination host.
        if !P.var_groups_high_functionality_express_earlier
            push!(dst_host.liver_infections, dst_inf)
        else
            genes_reorder = zeros(AlleleId, P.n_loci, P.n_genes_per_strain)
            group_ids = []
            for i in 1:size(dst_inf.genes)[2]
                gene_temp_alleles = dst_inf.genes[:,i]
                gene_temp = Gene(gene_temp_alleles)
                # @assert haskey(s.association_genes_to_var_groups, gene_temp)
                gene_temp_group_id = s.association_genes_to_var_groups[gene_temp]
                push!(group_ids, gene_temp_group_id)
            end
            
            index_temp_all = []
            for i in 1:length(func_rank)
                group_id_temp = findfirst(item->item == i, func_rank)
                index_temp = findall(item->item == group_id_temp, group_ids)
                push!(index_temp_all, length(index_temp))
                index_temp_reorder = if i == 1
                    1:length(index_temp)
                else 
                    (sum(index_temp_all[1:(i-1)])+1):sum(index_temp_all[1:i])
                end
                genes_reorder[:,index_temp_reorder] = dst_inf.genes[:, index_temp]
            end
            dst_inf.genes = genes_reorder
            
            push!(dst_host.liver_infections, dst_inf)
        end

        # Update population wide host liver max
        should_update_rates = should_update_rates || update_n_liver_infections_per_host_max(s, dst_host)
    end

    if should_update_rates
        # Only one event depends on n_liver_infections_per_host_max
        update_rate!(t, s, event_dist, LIVER_PROGRESS)
    end

    if transmitted
        stats.n_transmitting_bites += 1
        true
    else
        false
    end
end

"""
    Update maximum number of liver infections per host.

    Return whether the number was changed.
"""
function update_n_liver_infections_per_host_max(s, host)
    n_liver_infections = length(host.liver_infections)
    if s.n_liver_infections_per_host_max < n_liver_infections
        s.n_liver_infections_per_host_max = n_liver_infections
        true
    else
        false
    end
end

"""
    Update maximum number of active infections per host.

    Return whether the number was changed.
"""
function update_n_active_infections_per_host_max(s, host)
    n_active_infections = length(host.active_infections)
    if s.n_active_infections_per_host_max < n_active_infections
        s.n_active_infections_per_host_max = n_active_infections
        true
    else
        false
    end
end


function do_rebirth!(t, s, host)
    host.id = next_host_id!(s)
    host.t_birth = t
    # host.t_death = t + draw_host_lifetime()
    host.n_cleared_infections = 0
    empty!(host.liver_infections)
    empty!(host.active_infections)
    empty!(host.immunity)
end


### IMMIGRATION EVENT ###

function get_rate_immigration(t, s)
    P.immigration_rate_fraction * get_rate_biting(t, s) * s.infected_ratio
end

function do_immigration!(t, s, stats, event_dist)

    # Sample a random host and advance it (rebirth or infection activation).
    host = rand(s.rng, s.hosts)

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
    infection.liver_index = 1
    infection.expression_index = 0
    infection.expression_index_locus = 0
    if !P.var_groups_fix_ratio
        for i in 1:P.n_genes_per_strain
            infection.genes[:,i] = s.gene_pool[:, rand(s.rng, 1:size(s.gene_pool)[2])]
        end
        group_ids = []
        for i in 1:size(infection.genes)[2]
            group_id = s.association_genes_to_var_groups[Gene(infection.genes[:,i])]
            push!(group_ids, group_id)
        end
    else 
        infection_genes = zeros(AlleleId, P.n_loci, P.n_genes_per_strain)
        group_ids = []
        for group_id in 1:length(P.var_groups_ratio)
            if P.var_groups_ratio[group_id] > 0.0
                infection_genes_index_var_group = infection_genes_index_var_groups[group_id]
                for infection_gene_index_var_group in infection_genes_index_var_group
                    infection_gene = s.gene_pool[:, rand(s.rng, 1:size(s.gene_pool)[2])]
                    # @assert haskey(s.association_genes_to_var_groups, Gene(infection_gene))
                    infection_gene_group_id = s.association_genes_to_var_groups[Gene(infection_gene)]
                    while infection_gene_group_id != group_id
                        infection_gene = s.gene_pool[:, rand(s.rng, 1:size(s.gene_pool)[2])]
                        # @assert haskey(s.association_genes_to_var_groups, Gene(infection_gene))
                        infection_gene_group_id = s.association_genes_to_var_groups[Gene(infection_gene)]
                    end
                    # @assert s.association_genes_to_var_groups[Gene(infection_gene)] == group_id
                    push!(group_ids, group_id)
                    infection_genes[:,infection_gene_index_var_group] = infection_gene # sample with replacement
                end
            end
        end
        infection.genes[:,:] = reshape(infection_genes, (P.n_loci, P.n_genes_per_strain))
    end

    # Check if need to reorder genes when setting high functionality express earlier to be true.
    if P.var_groups_high_functionality_express_earlier
        genes_reorder = zeros(AlleleId, P.n_loci, P.n_genes_per_strain)
        index_temp_all = []
        for j in 1:length(func_rank)
            group_id_temp = findfirst(item->item == j, func_rank)
            index_temp = findall(item->item == group_id_temp, group_ids)
            push!(index_temp_all, length(index_temp))
            if length(index_temp) > 0
                index_temp_reorder = if j == 1
                    1:length(index_temp)
                else 
                    (sum(index_temp_all[1:(j-1)])+1):sum(index_temp_all[1:j])
                end
                genes_reorder[:,index_temp_reorder] = infection.genes[:, index_temp]
            end
        end
        infection.genes = genes_reorder    
    end

    # Add infection to host.
    push!(host.liver_infections, infection)
    
    if update_n_liver_infections_per_host_max(s, host)
        update_rate!(t, s, event_dist, LIVER_PROGRESS)
    end

    true
end


### RANDOM BACKGROUND CLEARANCE EVENT ###

function get_rate_background_clearance(t, s)
    P.background_clearance_rate * P.n_hosts * s.n_active_infections_per_host_max
end

function do_background_clearance(t, s, stats, event_dist)
    host = rand(s.rng, P.n_hosts)
    inf_index = rand(s.rng, 1:s.n_active_infections_per_host_max)

    # If the infection index is out of range, this is a rejected sample.
    # Otherwise we'll proceed.
    if inf_index > length(host.active_infections)
        false
    else
        infection = host.active_infections[inf_index]
        #println("do_background_clearance actually happening")
        clear_active_infection!(t, s, host, inf_index)
        true
    end
end


### LIVER PROGRESS EVENT ###

function get_rate_liver_progress(t, s)
    P.n_hosts * s.n_liver_infections_per_host_max * P.liver_erlang_shape / P.t_liver_stage
end

function do_liver_progress!(t, s, stats, event_dist)
    host = rand(s.rng, s.hosts)
    inf_index = rand(s.rng, 1:s.n_liver_infections_per_host_max)

    if inf_index > length(host.liver_infections)
        return false
    end
    infection = host.liver_infections[inf_index]

    # Increment liver progress or activate infection
    @assert infection.liver_index >= 1 && infection.liver_index <= P.liver_erlang_shape
    if infection.liver_index == P.liver_erlang_shape
        # If the infection is past the liver stage, remove it from the liver and activate it.
        delete_and_swap_with_end!(host.liver_infections, inf_index)
        # If there's room, move it into the active infections array.
        # Otherwise, just put it into the recycle bin.
        if isnothing(P.n_infections_active_max) || length(host.active_infections) < P.n_infections_active_max
            infection.liver_index = 0
            infection.expression_index = 1
            if P.whole_gene_immune
                infection.expression_index_locus = P.n_loci
            else
                infection.expression_index_locus = 1
            end
            push!(host.active_infections, infection)
            infection.t_expression = t
            (should_update_immunity_loss_rate, expression_ended) = advance_immune_genes!(t, s, host, length(host.active_infections))
            if should_update_immunity_loss_rate
                update_rate!(t, s, event_dist, IMMUNITY_LOSS)
            end

            # Update maximum number of active infections per host
            if update_n_active_infections_per_host_max(s, host)
                update_rates_from_n_active_infections_per_host_max!(t, s, event_dist)
            end
        else
            push!(s.old_infections, infection)
        end
    else
        infection.liver_index += 1
    end

    true
end


### SWITCHING EVENT ###

function get_rate_switching(t, s)
    # The total rate includes both active and liver infections because host state may not be fully up to date,
    # and a liver infection may be activated when host state is updated to the current time.
    # Rejection sampling is used to effect the correct rate.
    if !P.whole_gene_immune
        # Switching rate set by total number of alleles.
        (maximum(P.switching_rate) * P.n_loci) * P.n_hosts * s.n_active_infections_per_host_max
    else
        maximum(P.switching_rate) * P.n_hosts * s.n_active_infections_per_host_max
    end
end

function do_switching!(t, s, stats, event_dist)
    host = rand(s.rng, s.hosts)
    inf_index = rand(s.rng, 1:s.n_active_infections_per_host_max)
    
    # If the infection index is out of range, this is a rejected sample.
    # Otherwise we'll proceed.
    if inf_index > length(host.active_infections)
        return false
    end
    infection = host.active_infections[inf_index]

    gene_expression = infection.genes[:, infection.expression_index]
    # @assert haskey(s.association_genes_to_var_groups, Gene(gene_expression))
    gene_expression_group_id = s.association_genes_to_var_groups[Gene(gene_expression)]
    p_acceptance = P.switching_rate[gene_expression_group_id]/maximum(P.switching_rate)
    if rand(s.rng, ) < p_acceptance
        should_update_rates = false

        """
        Check that the host is not immune to the currently expressed loci of the gene
        """
        # immunity_level_current_gene_locus = get(host.immunity.vd[infection.expression_index_locus], infection.genes[infection.expression_index_locus, infection.expression_index], ImmunityLevel(0))
        # if immunity_level_current_gene_locus > 0
        #     println(immunity_level_current_gene_locus) 
        # end
        # @assert immunity_level_current_gene_locus == 0


        """
        Increment immunity level to currently expressed gene.
        For the partial allele model, expression advance by alleles, but
        Immunity only gains after the full gene finishes expression
        """
        if infection.expression_index_locus == P.n_loci
            should_update_rates = increment_immunity!(t, s, host, infection.genes[:, infection.expression_index])
        end

        # If we're at the end, clear the infection and return.
        if infection.expression_index == P.n_genes_per_strain && infection.expression_index_locus == P.n_loci
            clear_active_infection!(t, s, host, inf_index)
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
            (should_update_rates_i, expression_ended) = advance_immune_genes!(t, s, host, i)
            should_update_rates = should_update_rates || should_update_rates_i
            if !expression_ended
                # If there is no end of expression and reordering of infections, then index plus 1.
                i += 1
            end
        end

        if should_update_rates
            update_rate!(t, s, event_dist, IMMUNITY_LOSS)
        end

        true
    else
        false
    end
end

# This function moves the expression index of an infection to its first non-immune allele/gene.
function advance_immune_genes!(t, s, host, inf_index)
    # Advance expression until a non-immune gene or allele is reached.
    infection = host.active_infections[inf_index]
    # If the host not immune, stop advancing.
    should_advance = is_immune(host.immunity, infection.genes[:, infection.expression_index],infection.expression_index_locus)

    should_update_rates = false
    while should_advance
        # Increment immunity level to currently expressed gene or allele.
        if infection.expression_index_locus == P.n_loci
            should_update_rates = should_update_rates || increment_immunity!(t, s, host, infection.genes[:, infection.expression_index])
        end

        # If we're at the end, clear the infection and return.
        if infection.expression_index == P.n_genes_per_strain && infection.expression_index_locus == P.n_loci
            clear_active_infection!(t, s, host, inf_index)

            # Return whether rates need to be updated, and `true` to indicate that expression ended
            return (should_update_rates, true)
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
        should_advance = is_immune(host.immunity, infection.genes[:, infection.expression_index],infection.expression_index_locus)
    end
    # Return whether to update immunity loss rate, and `false` to indicate that expression did not end
    (should_update_rates, false)
end


### MUTATION EVENT ###
# Update mutation and recombination rates towards all infections.
function get_rate_mutation(t, s)
    # The total rate includes both active and liver infections because host state may not be fully up to date,
    # and a liver infection may be activated when host state is updated to the current time.
    # Rejection sampling is used to effect the correct rate.
    P.mutation_rate * P.n_hosts * s.n_active_infections_per_host_max * P.n_genes_per_strain * P.n_loci
end

function do_mutation!(t, s, stats, event_dist)
    host = rand(s.rng, s.hosts)
    inf_index = rand(s.rng, 1:s.n_active_infections_per_host_max)

    # If there's no active infection at the drawn index, reject this sample.
    if inf_index > length(host.active_infections)
        return false
    end

    expression_index = rand(s.rng, 1:P.n_genes_per_strain)
    locus = rand(s.rng, 1:P.n_loci)

    infection = host.active_infections[inf_index]

    # If we ever generate too many alleles for 16-bit ints, we'll need to use bigger ones.
    @assert s.n_alleles[locus] < typemax(AlleleId)

    # Generate a new allele and insert it at the drawn location.
    s.n_alleles[locus] += 1
    # Add new gene and its group id to the association map. Mutation preserves group id. 
    gene_ori = Gene(infection.genes[:, expression_index])
    infection.genes[locus, expression_index] = s.n_alleles[locus]
    gene_mut = Gene(infection.genes[:, expression_index])
    # @assert haskey(s.association_genes_to_var_groups, gene_ori) # Check the association map contains the group id of the original gene, which should always be true.
    gene_mut_group_id = s.association_genes_to_var_groups[gene_ori]
    # @assert !(haskey(s.association_genes_to_var_groups, gene_mut))
    s.association_genes_to_var_groups[gene_mut] = gene_mut_group_id
    # @assert haskey(s.association_genes_to_var_groups, gene_mut) # Make sure the group id for the mutated gene is added to the association map.
    infection.strain_id = next_strain_id!(s)
    s.next_gene_id_mut += 1

    true
end


### ECTOPIC RECOMBINATION EVENT ###

function get_rate_ectopic_recombination(t, s)
    # The total rate includes both active and liver infections because host state may not be fully up to date,
    # and a liver infection may be activated when host state is updated to the current time.
    # Rejection sampling is used to effect the correct rate.
    maximum(P.ectopic_recombination_rate)^2 *
        P.n_hosts * s.n_active_infections_per_host_max *
        P.n_genes_per_strain * (P.n_genes_per_strain - 1) / 2.0
end



function do_ectopic_recombination!(t, s, stats, event_dist)
    host = rand(s.rng, s.hosts)
    inf_index = rand(s.rng, 1:s.n_active_infections_per_host_max)

    # If there's no active infection at the drawn index, reject this sample.
    if inf_index > length(host.active_infections)
        return false
    end

    infection = host.active_infections[inf_index]

    gene_indices = rand(s.rng, gene_pair_indices)
    gene_index_1 = gene_indices[1]
    gene_index_2 = gene_indices[2]

    gene1 = infection.genes[:, gene_index_1]
    gene2 = infection.genes[:, gene_index_2]

    # Check the group id of the two drawn genes, reject with a probability accounting for the differential ectopic recombination rate across var groups
    gene1_group_id = s.association_genes_to_var_groups[Gene(gene1)]
    gene2_group_id = s.association_genes_to_var_groups[Gene(gene2)]
    
    # if different var groups do not share alleles, then genes from different groups do not recombine either. 
    if P.var_groups_do_not_share_alleles && gene1_group_id != gene2_group_id
        return false
    end
    
    gene1_recomb_rate_ratio = P.ectopic_recombination_rate[gene1_group_id]/maximum(P.ectopic_recombination_rate)
    gene2_recomb_rate_ratio = P.ectopic_recombination_rate[gene2_group_id]/maximum(P.ectopic_recombination_rate)
    p_acceptance = gene1_recomb_rate_ratio * gene2_recomb_rate_ratio

    if rand(s.rng) < p_acceptance
        # If the genes are the same, this is a no-op.
        if gene1 == gene2
            return false
        end

        breakpoint, p_functional = if P.ectopic_recombination_generates_new_alleles
            # Choose a breakpoint.
            breakpoint = P.n_loci * rand(s.rng)
            p_functional = p_recombination_is_functional_real(gene1, gene2, breakpoint)

            (Int(ceil(breakpoint)), p_functional)
        else
            # Choose a breakpoint.
            breakpoint = rand(s.rng, 1:P.n_loci)
            if breakpoint == 1
                return false
            end

            p_functional = p_recombination_is_functional_integer(gene1, gene2, breakpoint)

            (breakpoint, p_functional)
        end

        is_conversion = rand(s.rng) < P.p_ectopic_recombination_is_conversion

        recombined = false

        create_new_allele = P.ectopic_recombination_generates_new_alleles &&
            rand(s.rng) < P.p_ectopic_recombination_generates_new_allele

        # Recombine to modify first gene, if functional.
        if !is_conversion && rand(s.rng) < p_functional
            infection.genes[:, gene_index_1] = if create_new_allele
                recombine_genes_new_allele(s, gene1, gene2, breakpoint)
            else
                recombine_genes(gene1, gene2, breakpoint, s)
            end
        end

        # Recombine to modify second gene, if functional.
        if rand(s.rng) < p_functional
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
            if gene1 != infection.genes[:, gene_index_1]
                gene1_ectopic = Gene(infection.genes[:, gene_index_1])
                # @assert haskey(s.association_genes_to_var_groups, Gene(gene1))
                gene1_ectopic_group_id = s.association_genes_to_var_groups[Gene(gene1)]
                # @assert !(haskey(s.association_genes_to_var_groups, gene1_ectopic))
                if !(haskey(s.association_genes_to_var_groups, gene1_ectopic))
                    s.association_genes_to_var_groups[gene1_ectopic] = gene1_ectopic_group_id
                end
                # @assert haskey(s.association_genes_to_var_groups, gene1_ectopic)
            end
            if gene2 != infection.genes[:, gene_index_2]
                gene2_ectopic = Gene(infection.genes[:, gene_index_2])
                # @assert haskey(s.association_genes_to_var_groups, Gene(gene2))
                gene2_ectopic_group_id = s.association_genes_to_var_groups[Gene(gene2)]
                # @assert !(haskey(s.association_genes_to_var_groups, gene2_ectopic))
                if !(haskey(s.association_genes_to_var_groups, gene2_ectopic))
                    s.association_genes_to_var_groups[gene2_ectopic] = gene2_ectopic_group_id
                end
                # @assert haskey(s.association_genes_to_var_groups, gene2_ectopic)
            end
            infection.strain_id = next_strain_id!(s)
            true
        else
            false
        end
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
    if P.sample_infection_duration_every !== nothing
        if s.n_cleared_infections % P.sample_infection_duration_every == 0
            # Calculate and write the infection duration.
            add_infection_duration!(t, s, host, inf_index)
        end
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

function do_immunity_loss!(t, s, stats, event_dist)
    host = rand(s.rng, s.hosts)
    immunity_index = rand(s.rng, 1:s.n_immunities_per_host_max)

    # If the immunity index is beyond this host's immunity count, reject this sample.
    if immunity_index >  immunity_count(host.immunity)
        false
    else
        decrement_immunity_at_sampled_index!(host.immunity, immunity_index)
        true
    end
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

function empty!(ih::ImmuneHistory)
    for d in ih.vd
        empty!(d)
    end
end

function immunity_count(ih::ImmuneHistory)
    sum(length(ih.vd[locus]) for locus in 1:length(ih.vd))
end


function increment_immunity!(t, s, host, gene)
    increment_immunity!(host.immunity, gene)

    n_immunities = immunity_count(host.immunity)
    if s.n_immunities_per_host_max < n_immunities
        s.n_immunities_per_host_max = n_immunities
        true
    else
        false
    end
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
