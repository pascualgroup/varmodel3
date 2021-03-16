function run_discrete_time(p::Params)
    s = State(p)
    db = initialize_database(p)

    start_datetime = now()
    last_summary_datetime = start_datetime
    n_infected_bites = 0
    n_infected_bites_with_space = 0
    n_total_bites = 0

    execute(db, "BEGIN TRANSACTION")
    write_summary(0, p, s, db, missing, missing, missing, missing)
    execute(db, "COMMIT")

    verify(p, s)
    
    execute(db, "BEGIN TRANSACTION")

    # TODO: use dt
    for t in 1:p.t_end
#         println("stepping to t = $(t)")

        do_rebirth!(t, p, s)
        do_immunity_loss!(t, p, s)
        do_activation!(t, p, s)
        do_switching!(t, p, s)

        n_infected_bites_t, n_infected_bites_with_space_t, n_total_bites_t = do_biting!(t, p, s)
        n_infected_bites += n_infected_bites_t
        n_infected_bites_with_space += n_infected_bites_with_space_t
        n_total_bites += n_total_bites_t

        do_mutation!(t, p, s)
        do_recombination!(t, p, s)
        do_immigration!(t, p, s)

        if t % p.summary_period == 0
            execute(db, "COMMIT")

            execute(db, "BEGIN TRANSACTION")

            next_summary_datetime = now()
            write_summary(
                t, p, s, db,
                Dates.value(next_summary_datetime - last_summary_datetime),
                n_infected_bites, n_infected_bites_with_space, n_total_bites
            )
            last_summary_datetime = next_summary_datetime
            n_infected_bites = 0
            n_infected_bites_with_space = 0
            n_total_bites = 0
        end

        if t % p.strain_count_period == 0
            write_gene_strain_counts(t, p, s, db)
        end

        if t % p.summary_period == 0
            execute(db, "COMMIT")
            execute(db, "BEGIN TRANSACTION")
        end

        if p.verification_period != nothing && t % p.verification_period == 0
            verify(p, s)
        end
    end
end

function write_summary(t, p, s, db, elapsed_time_ms, n_infected_bites, n_infected_bites_with_space, n_total_bites)
    infections_liver = .!isnan.(s.t_infection_liver)
    n_infections_liver = sum(infections_liver)

    infections_active = s.expression_index .!= 0
    n_infections_active = sum(infections_active)

    n_infections = n_infections_liver + n_infections_active

    infected_liver = sum(infections_liver, dims = 2) .!= 0
    infected_active = sum(infections_active, dims = 2) .!= 0
    n_infected_liver = sum(infected_liver)
    n_infected_active = sum(infected_active)
    n_infected = sum(infected_liver .| infected_active)

    exec_time = elapsed_time_ms

    println("")
    println("t = $(t)")
    println("n_infections = $(n_infections)")

    execute(db.summary, [
        t,
        n_infections_liver,
        n_infections_active,
        n_infections,
        n_infected_liver,
        n_infected_active,
        n_infected,
        n_infected_bites,
        n_infected_bites_with_space,
        n_total_bites,
        exec_time
    ])
end

function write_gene_strain_counts(t, p, s, db)
    genes::Set{SVector{p.n_loci, AlleleId}} = Set()
    strains::BitSet = BitSet()
    
    count_circulating_genes_and_strains!(p, s.strain_id_liver, s.genes_liver, genes, strains)
    count_circulating_genes_and_strains!(p, s.strain_id_active, s.genes_active, genes, strains)
    
    execute(db.gene_strain_counts, [t, length(genes), length(strains)])
end

function count_circulating_genes_and_strains!(p, infection_strain_id, infection_genes, genes, strains)
    for k in 1:size(infection_strain_id)[2]
        for j in 1:size(infection_strain_id)[1]
            strain_id = infection_strain_id[j, k]
            if strain_id != 0 && !(strain_id in strains)
                for i in 1:size(infection_genes)[2]
                    push!(genes, @view(infection_genes[:, i, j, k]))
                end
                push!(strains, strain_id)
            end
        end
    end
end

function do_rebirth!(t, p, s)
#     println("do_rebirth!()")

    dead_hosts = findall(s.t_death .< t)
#     println("n dead: $(length(dead_hosts))")

    s.t_birth[dead_hosts] = s.t_death[dead_hosts]
    s.t_death[dead_hosts] = s.t_birth[dead_hosts] + [draw_host_lifetime(p) for i in 1:length(dead_hosts)]
    s.t_infection_liver[:, dead_hosts] .= NaN32
    s.t_infection_active[:, dead_hosts] .= NaN32
    s.strain_id_liver[:, dead_hosts] .= 0
    s.strain_id_active[:, dead_hosts] .= 0
    s.genes_liver[:, :, :, dead_hosts] .= 0
    s.genes_active[:, :, :, dead_hosts] .= 0
    s.expression_index[:, dead_hosts] .= 0
    
    for host in dead_hosts
        empty!(s.immunity[host])
    end
    
#     verify(p, s)
end

function do_immunity_loss!(t, p, s)
#     println("do_immunity_loss!()")
    
    p_loss = 1 - exp(-p.immunity_loss_rate)
    n_immunity_max = maximum(length(s.immunity[i]) for i in 1:p.n_hosts)
    n_loss = rand(Binomial(n_immunity_max * p.n_hosts))
    
    indices = CartesianIndices((n_immunity_max, p.n_hosts))[
        sample(1:(n_immunity_max * p.n_hosts), n_loss)
    ]
#     println("n_loss_raw = $(length(indices))")
    
    for index in indices
        immunity_index = index[1]
        host = index[2]
        
        # If the immunity index is greater than the number of immune genes,
        # do nothing (rejection sampling)
        # Otherwise, select a key via iteration order (inefficient but maybe not bottleneck...)
        if immunity_index <= length(s.immunity[host])
#             println("losing $(host), $(index)")
            adjust_immunity!(p, s, host, get_key_by_iteration_order(s.immunity[host], immunity_index), -1)
        end
    end
    
#     verify(p, s)
end

function findfirst_stable(x)
    y = findfirst(x)
    y == nothing ? 0 : y
end

function do_activation!(t, p, s)
#     println("do_activation!()")
    
    # Identify hosts with infections ready to become active
    ready_mask_raw = s.t_infection_liver .+ p.t_liver_stage .<= t
    host_indices = findall(
        reshape(sum(ready_mask_raw, dims = 1), :) .> 0
    )
    ready_mask = ready_mask_raw[:, host_indices]
    ready_count = reshape(sum(ready_mask, dims = 1), :)

    # Identify available slots
    available_mask = s.expression_index[:, host_indices,] .== 0
    available_count = reshape(sum(available_mask, dims = 1), :)

    # Compute the number that will actually be activated
    activation_count = min.(ready_count, available_count)

    # Match indices from liver stage arrays to active arrays
    src_mask = mask_with_column_limits(ready_mask, activation_count)
    src_indices = findall(src_mask)
    dst_mask = mask_with_column_limits(available_mask, activation_count)
    dst_indices = findall(dst_mask)

    @assert length(src_indices) == length(dst_indices)

    expression_index = @view(s.expression_index[:, host_indices])
    t_infection_liver = @view(s.t_infection_liver[:, host_indices])
    t_infection_active = @view(s.t_infection_active[:, host_indices])
    strain_id_liver = @view s.strain_id_liver[:, host_indices]
    strain_id_active = @view s.strain_id_active[:, host_indices]
    genes_liver = @view s.genes_liver[:, :, :, host_indices]
    genes_active = @view s.genes_active[:, :, :, host_indices]

    # Update active infection data
    expression_index[dst_indices] .= 1
    t_infection_active[dst_indices] = t_infection_liver[src_indices]

    strain_id_active[dst_indices] = strain_id_liver[src_indices]
    genes_active[:, :, dst_indices] = genes_liver[:, :, src_indices]

    # Deactivate all ready infections (including those that didn't get activated)
    t_infection_liver[src_indices] .= NaN32
    strain_id_liver[src_indices] .= 0
    genes_liver[:, :, src_indices] .= 0

#     verify(p, s)
end

function do_switching!(t, p, s)
#     println("do_switching!()")
    # TODO: adjust for dt
    p_switch = 1 - exp(-p.switching_rate)
    n_infections = p.n_hosts * p.n_infections_active_max

    n_trans_raw = rand(Binomial(n_infections, p_switch))
    if n_trans_raw == 0
        return
    end

    # Uniformly randomly sample indices in one-dimensional array
    # Includes nonexistent infections
    indices_raw = sample(1:n_infections, n_trans_raw, replace = false)
    
    # Filter indices down to active infections
    indices = CartesianIndices(s.expression_index)[indices_raw[s.expression_index[indices_raw] .> 0]]
    if length(indices) == 0
        return
    end
    
    # Loop in parallel through advancing expression for all sampled infections
    # until they have all reached a non-immune gene
    expression_index = @view s.expression_index[indices]
    t_infection_active = @view s.t_infection_active[indices]
    strain_id_active = @view s.strain_id_active[indices]
    genes_active = @view s.genes_active[:,:,indices]
    
    done = falses(length(indices))
    while true
        # Increment immunity level for all infections
        for index in indices[.!done]
            inf = index[1]
            host = index[2]
            adjust_immunity!(p, s, host, @view(s.genes_active[:, s.expression_index[index], inf, host]), 1)
        end
        
        # Clear infections at last index
        to_clear = expression_index .== p.n_genes_per_strain
        t_infection_active[to_clear] .= NaN32
        strain_id_active[to_clear] .= 0
        genes_active[:, :, to_clear] .= 0
        expression_index[to_clear] .= 0
        done[to_clear] .= true
        
        # Advance expression for infections not yet done
        not_done = .!done
        expression_index[not_done] .+= 1
        
        # Stop advancing if not immune
        for (i, index) in enumerate(indices)
            if !done[i]
                inf = index[1]
                host = index[2]
                if !is_immune(p, s, host, (s.genes_active[:, s.expression_index[index], inf, host]))
                    done[i] = true
                end
            end
        end
        
        if all(done)
            break
        end
    end
    
#     verify(p, s)
end

function is_immune(p, s, host, gene)
    get(s.immunity[host], gene, 0) > 0
end

function adjust_immunity!(p, s, host, gene, amount)
    @assert amount != 0
    
    # Get old and new immunity level
    old_level = get(s.immunity[host], gene, ImmunityLevel(0))
    new_level = clamp(old_level + amount, 0, p.immunity_level_max)
    
    # Update dictionaries
    if new_level == 0
        delete!(s.immunity[host], gene)
    else
        s.immunity[host][gene] = new_level
    end
    nothing
end

function do_biting!(t, p, s)
#     println("do_biting!()")

    # TODO: adjust for dt
    p_bite = 1 - exp(-p.biting_rate[1 + Int(floor(t)) % p.t_year])
    n_bites_raw = rand(Binomial(p.n_hosts, p_bite))

    # Sample source hosts and destination hosts
    src_hosts_raw = sample(1:p.n_hosts, n_bites_raw, replace = false)
    dst_hosts_raw = sample(1:p.n_hosts, n_bites_raw, replace = false)
    if length(src_hosts_raw) == 0
        return (0, 0, 0)
    end

    # Identify indices with active source infections and available space in destination
    src_valid = reshape(any(s.expression_index[:,src_hosts_raw] .> 0, dims = 1), n_bites_raw)
    dst_valid = reshape(any(isnan.(s.t_infection_liver[:,dst_hosts_raw]), dims = 1), n_bites_raw)
    src_dst_valid = src_valid .& dst_valid
    
    # Record number of infected bites for output
    n_infected_bites = sum(src_valid)
    
    # Filter source and destination host indices
    src_hosts = src_hosts_raw[src_dst_valid]
    if length(src_hosts) == 0
        return (0, 0, 0)
    end
    dst_hosts = dst_hosts_raw[src_dst_valid]
    
    # Record number of infected bites with space in destination for output
    n_bites = length(dst_hosts)
#     println("n_bites = $(n_bites)")

    # Compute number of active source infections and their locations in array
    src_active_mask = s.expression_index[:, src_hosts] .> 0
    src_inf_count = sum(src_active_mask; dims = 1) # dims: (1, n_bites)
    
    # Compute probability of infection, dims: (1, n_bites)
    p_infect = if p.coinfection_reduces_transmission
        fill(p.transmissibility, (1, n_bites))
    else
        p.transmissibility ./ src_inf_count
    end
    
    # Identify which source infections will be involved in transmission
    src_transmit_mask = src_active_mask .& (rand(Float32, size(src_active_mask)) .< p_infect)
    src_transmit_count = reshape(sum(src_transmit_mask; dims = 1), :)
#     println("src_transmit_count = $(src_transmit_count)")
    
    # Identify available slots in destination livers
    dst_avail_mask = s.strain_id_liver[:, dst_hosts] .== 0
    dst_avail_count = reshape(sum(dst_avail_mask; dims = 1), :)
#     println("dst_avail_count = $(dst_avail_count)")

    # Compute how many infections will be transmitted in each bite
    transmit_count = min.(src_transmit_count, dst_avail_count)
#     println("transmit_count = $(transmit_count)")
    transmit_count_total = sum(transmit_count)
    
    # Do one transmission at a time in parallel across hosts
    for i in 1:maximum(transmit_count)
        bite_mask = (transmit_count .> 0) .& (i .<= transmit_count)
        n_bites_i = sum(bite_mask)
        
        src_hosts_i = src_hosts[bite_mask]
        dst_hosts_i = dst_hosts[bite_mask]
        
        # two vectors of host-specific indices for infections being recombined
        # dimensions: (n_bites_i,)
        src_transmit_indices_1 = sample_true_indices_by_column(src_transmit_mask[:, bite_mask])
        src_transmit_indices_2 = sample_true_indices_by_column(src_transmit_mask[:, bite_mask])
    
        # two vectors of strain IDs for strains being recombined
        # dimensions: (n_bites_i,)
        strain_ids_1 = s.strain_id_active[zip_cartesian(src_transmit_indices_1, src_hosts_i)]
        strain_ids_2 = s.strain_id_active[zip_cartesian(src_transmit_indices_2, src_hosts_i)]
        
        @assert all(strain_ids_1 .> 0)
        @assert all(strain_ids_2 .> 0)
        
        # Genes being transmitted
        genes_1 = s.genes_active[:, :, zip_cartesian(src_transmit_indices_1, src_hosts_i)]
        genes_2 = s.genes_active[:, :, zip_cartesian(src_transmit_indices_2, src_hosts_i)]
    
        # Recombined strain IDs.
        # Reuse IDs for paired strains with the same ID. Make new IDs for the other pairs.
        # Dimensions: (n_bites_i,)
        same_strain = strain_ids_1 .== strain_ids_2
        not_same_strain = .!same_strain
        n_new = sum(not_same_strain)
        strain_ids_recombined =
            (same_strain .* strain_ids_1) .+
            fill_mask_with_entries(not_same_strain, next_strain_ids!(s, n_new))
    
        # Recombined genes.
        # For non-recombining infections (given by same_strain), shuffle by column.
        # For recombining infections, sample from genes_1 and genes_2 stacked on top of each other.
        genes_recombined = fill(0, (p.n_loci, p.n_genes_per_strain, n_bites_i))
        genes_recombined[:, :, same_strain] = shuffle_genes_within_strains(
            genes_1[:, :, same_strain]
        )
        genes_recombined[:, :, not_same_strain] = recombine_multiple_strains(
            genes_1[:, :, not_same_strain],
            genes_2[:, :, not_same_strain]
        )
        
        # Copy in new infections
        dst_transmit_indices = findfirst_each_column(s.strain_id_liver[:,dst_hosts_i] .== 0)
        inf_host_indices = zip_cartesian(dst_transmit_indices, dst_hosts_i)
        s.t_infection_liver[inf_host_indices] .= t
        s.strain_id_liver[inf_host_indices] = strain_ids_recombined
        s.genes_liver[:, :, inf_host_indices] = genes_recombined
    end
    
#     verify(p, s)
    
    (n_infected_bites, n_bites, n_bites_raw)
end

function shuffle_genes_within_strains(genes)
    @assert length(size(genes)) == 3
    
    genes_shuffled = Array(genes)
    for i in 1:size(genes)[3]
        genes_shuffled[:,:,i] = genes_shuffled[:,randperm(size(genes)[2]),i]
    end
    genes_shuffled
end

function recombine_multiple_strains(genes_1, genes_2)
    all_genes = cat(genes_1, genes_2; dims = 2)
    genes_recombined = similar(genes_1)
    n_genes = size(genes_1)[2]
    for i in 1:size(genes_1)[3]
        genes_recombined[:,:,i] = all_genes[:, sample(1:(n_genes * 2), n_genes, replace = false), i]
    end
    genes_recombined
end

function do_mutation!(t, p, s)
    # TODO: adjust for dt
    p_mut = 1 - exp(-p.mutation_rate)

    # Draw mutations randomly using linear indexing
    n_possible_mut = length(s.genes_active)
    n_mut_raw = rand(Binomial(n_possible_mut, p_mut))
    mut_indices_raw = sample(1:n_possible_mut, n_mut_raw, replace = false)

    # Apply mutations for infections that actually exist
    genes_cartesian = CartesianIndices(s.genes_active)
    mut_indices = filter(
        function(i)
            ind = genes_cartesian[i] # Gets multidimensional index from linear index
            inf = ind[3]
            host = ind[4]
            s.expression_index[inf, host] > 0
        end,
        mut_indices_raw
    )
    n_mut = length(mut_indices)
#     println("n_mut_raw = $(n_mut_raw), n_mut = $(n_mut)")

    # Apply mutations by creating new alleles
    for i in mut_indices
        (locus, exp_ind, inf, host) = Tuple(genes_cartesian[i])
        @assert s.expression_index[inf, host] > 0
        @assert maximum(s.n_alleles) < typemax(AlleleId)
        s.n_alleles[locus] += 1
        s.genes_active[locus, exp_ind, inf, host] = s.n_alleles[locus]

        s.strain_id_active[inf, host] = next_strain_id!(s)
    end
    
#     verify(p, s)
end

function do_recombination!(t, p, s)
#     println("do_recombination!()")

    p_recomb = 1 - exp(
        -p.ectopic_recombination_rate * (p.n_genes_per_strain) * (p.n_genes_per_strain - 1) / 2.0
    )

    n_inf_raw = p.n_hosts * p.n_infections_active_max
    n_recomb_raw = rand(Binomial(n_inf_raw, p_recomb))
    recomb_indices_raw = sample(1:n_inf_raw, n_recomb_raw, replace = false)

    # Identify active infections
    recomb_indices = filter(
        i -> s.expression_index[i] > 0,
        recomb_indices_raw
    )
    n_recomb = length(recomb_indices)

#     println("n_recomb_raw = $(n_recomb_raw), n_recomb = $(n_recomb)")

    # Apply recombinations
    # TODO: data-parallelize?
    inf_cartesian = CartesianIndices(s.expression_index)
    breakpoints = rand(1:p.n_loci, n_recomb)
    for i in 1:n_recomb
        (inf, host) = Tuple(inf_cartesian[recomb_indices[i]])

        i1, i2 = samplepair(p.n_genes_per_strain)
        src_gene_1 = @view s.genes_active[:, i1, inf, host]
        src_gene_2 = @view s.genes_active[:, i2, inf, host]
        if src_gene_1 == src_gene_2
            continue
        end

        breakpoint = breakpoints[i]
        if breakpoint == 1
            # Do nothing
        else
            similarity = gene_similarity(p, src_gene_1, src_gene_2, breakpoint)

            new_gene_1 = if rand() < similarity
                recombine_genes(src_gene_1, src_gene_2, breakpoint)
            else
                src_gene_1
            end

            new_gene_2 = if rand() < similarity
                recombine_genes(src_gene_2, src_gene_1, breakpoint)
            else
                src_gene_2
            end

            src_gene_1[:] = new_gene_1
            src_gene_2[:] = new_gene_2

            s.strain_id_active[inf, host] = s.next_strain_id
            s.next_strain_id += 1
        end
    end
    
#     verify(p, s)
end

function do_immigration!(t, p, s)
#     println("do_immigration!()")

    # TODO: adjust for dt
    p_bite = 1 - exp(
        -p.immigration_rate_fraction * p.biting_rate[1 + Int(floor(t)) % p.t_year]
    )
#     println("p_bite = $(p_bite)")
    n_bites_raw = rand(Binomial(p.n_hosts, p_bite))
    
    # Sample source hosts with available infection slots
    hosts_raw = sample(1:p.n_hosts, n_bites_raw, replace = false)
    valid = reshape(any(s.strain_id_liver[:,hosts_raw] .== 0, dims = 1), n_bites_raw)
    n_bites = sum(valid)
    
    if n_bites == 0
        return
    end
    
    hosts = hosts_raw[valid]
    
    inf_inds = findfirst_each_column(s.strain_id_liver[:,hosts] .== 0)
    inf_host_inds = zip_cartesian(inf_inds, hosts)
    
    s.t_infection_liver[inf_host_inds] .= t
    s.strain_id_liver[inf_host_inds] = next_strain_ids!(s, n_bites)
    genes = reshape(
        s.gene_pool[:, rand(1:size(s.gene_pool)[2], n_bites * p.n_genes_per_strain)],
        (p.n_loci, p.n_genes_per_strain, n_bites)
    )
    s.genes_liver[:, :, inf_host_inds] = genes
    
#     verify(p, s)
end
