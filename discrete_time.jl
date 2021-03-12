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
    println("n_infections = $(sum(s.expression_index .> 0))")

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

function count_circulating_genes_and_strains!(p::Params, infection_strain_id, infection_genes, genes, strains)
    for j in 1:size(infection_strain_id)[2]
        for i in 1:size(infection_strain_id)[1]
            strain_id = infection_strain_id[i,j]
            if strain_id != 0 && !(strain_id in strains)
                for k in 1:size(infection_genes)[3]
                    push!(genes, @view(infection_genes[i, j, k, :]))
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
    s.genes_liver[:, :, dead_hosts] .= 0
    s.genes_active[:, :, dead_hosts] .= 0
    s.expression_index[:, dead_hosts] .= 0
    # TODO: reset immunity
    
    verify(p, s)
end

function do_immunity_loss!(t, p, s)
#     println("do_immunity_loss!()")
    # TODO
end

function findfirst_stable(x)
    y = findfirst(x)
    y == nothing ? 0 : y
end

function do_activation!(t, p, s)
#     println("do_activation!()")
    
    # Identify hosts with infections ready to become active
    ready_mask_raw = s.t_infection_liver .+ p.t_liver_stage .<= t
    host_indices = findall(reshape(sum(ready_mask_raw, dims = 2), :) .> 0)
    ready_mask = ready_mask_raw[host_indices, :]
    ready_count = reshape(sum(ready_mask, dims = 2), :)

    # Identify available slots
    available_mask = s.expression_index[host_indices, :] .== 0
    available_count = reshape(sum(available_mask, dims = 2), :)

    # Compute the number that will actually be activated
    activation_count = min.(ready_count, available_count)

    # Match indices from liver stage arrays to active arrays
    src_mask = mask_with_row_limits(ready_mask, activation_count)
    src_indices = [CartesianIndex(x[2], x[1]) for x in findall(transpose(src_mask))]
    dst_mask = mask_with_row_limits(available_mask, activation_count)
    dst_indices = [CartesianIndex(x[2], x[1]) for x in findall(transpose(dst_mask))]

    @assert length(src_indices) == length(dst_indices)

    expression_index = @view(s.expression_index[host_indices, :])
    t_infection_liver = @view(s.t_infection_liver[host_indices, :])
    t_infection_active = @view(s.t_infection_active[host_indices, :])
    strain_id_liver = @view s.strain_id_liver[host_indices, :]
    strain_id_active = @view s.strain_id_active[host_indices, :]
    genes_liver = @view s.genes_liver[host_indices, :, :, :]
    genes_active = @view s.genes_active[host_indices, :, :, :]

    # Update active infection data
    expression_index[dst_indices] .= 1
    t_infection_active[dst_indices] = t_infection_liver[src_indices]

    strain_id_active[dst_indices] = strain_id_liver[src_indices]
    genes_active[dst_indices, :, :] = genes_liver[src_indices, :, :]

    # Deactivate all ready infections (including those that didn't get activated)
    t_infection_liver[src_indices] .= NaN32
    strain_id_liver[src_indices] .= 0
    genes_liver[src_indices, :, :] .= 0

    verify(p, s)
end

function do_switching!(t, p, s)
#     println("do_switching!()")
    # TODO: adjust for dt
    p_switch = 1 - exp(-p.switching_rate)

    n_infections = p.n_hosts * p.n_infections_active_max

    n_trans_raw = rand(Binomial(n_infections, p_switch))
#     println("n_trans_raw: $(n_trans_raw)")
    if n_trans_raw == 0
        return
    end

    # Uniformly randomly sample indices in one-dimensional array
    # Includes nonexistent infections
    indices_raw = sample(1:n_infections, n_trans_raw, replace = false)

    # Get one-dimensional views on infection
    expression_index_view = reshape(s.expression_index, n_infections)
    infection_genes_view = reshape(s.genes_active, (p.n_genes_per_strain, n_infections))

    # Filter indices down to active infections
    indices = indices_raw[expression_index_view[indices_raw] .> 0]
    n_trans = length(indices)
#     println("n_trans = $(n_trans)")
    if n_trans == 0
        return
    end

    # Compute host and infection indices from sampled linear indices
    cartesian_indices = cartesian_indices_infection_host(s)[indices]

    # Advance expression for each chosen infection
    for ci in cartesian_indices
        inf_index = ci[1]
        host_index = ci[2]

        # Advance expression until non-immune gene is reached
        while true
            exp_index = s.expression_index[ci]
            @assert exp_index > 0

            # TODO: gain extra immunity to current gene

            if exp_index == p.n_genes_per_strain
                # Clear infection if we're at the end
                s.t_infection_active[ci] = NaN32
                s.strain_id_active[ci] = 0
                s.genes_active[:, ci] .= 0
                s.expression_index[ci] = 0
                break
            else
                # Advance expression if we're not yet at the end
                s.expression_index[ci] += 1

                # TODO: If we're not immune to this gene, stop advancing expression
                # in the meantime: always stop advancing.
                break
            end
        end
    end
    
    verify(p, s)
end

# TODO: rewrite this whole function to be cleaner and simpler
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
    src_transmit_count = sum(src_transmit_mask; dims = 1)
    
    # Identify available slots in destination livers
    dst_avail_mask = s.strain_id_liver[:, dst_hosts] .> 0
    dst_avail_count = sum(dst_avail_mask; dims = 1) # dims: (1, n_bites)

    # Compute how many infections will be transmitted in each bite
    transmit_count = reshape(min.(src_transmit_count, dst_avail_count), n_bites)
    transmit_count_total = sum(transmit_count)
    
    # Do one transmission at a time in parallel across hosts
    for i in 1:maximum(transmit_count)
        bite_mask = transmit_count .<= i
        n_bites_i = sum(bite_mask)
        
        src_hosts_i = src_hosts[bite_mask]
        dst_hosts_i = dst_hosts[bite_mask]
        
        # two vectors of host-specific indices for infections being recombined
        # dimensions: (n_bites_i,)
        src_transmit_indices_1 = sample_true_indices_by_column(src_transmit_mask[:, bite_mask])
        src_transmit_indices_2 = sample_true_indices_by_column(src_transmit_mask[:, bite_mask])
    
        # two vectors of strain IDs for strains being recombined
        # dimensions: (n_bites_i,)
        strain_ids_1 = zip_index(s.strain_id_active, src_transmit_indices_1, src_hosts_i)
        strain_ids_2 = zip_index(s.strain_id_active, src_transmit_indices_2, src_hosts_i)
    
        @assert all(strain_ids_1 .> 0)
        @assert all(strain_ids_2 .> 0)
    
        # host-specific indices of genes for strains being recombined
        # first need indices of infections being transmitted
        # dimensions: (n_genes_per_strain, n_bites_i)
        gene_inds_1 = s.genes_active[:, zip_cartesian(src_transmit_indices_1, src_hosts_i)]
        gene_inds_2 = s.genes_active[:, zip_cartesian(src_transmit_indices_2, src_hosts_i)]
    
        # extract actual genes for strains being recombined
        # (1) treat gene_inds as long vectors of stacked columns to extract genes from host_genes
        # (2) reshape them into arrays corresponding to sequences of genes within infections
        # dimensions: (n_loci, n_genes_per_strain, n_bites_i)
        genes_1 = reshape(s.host_genes[:, gene_inds_1[:]], (p.n_loci, p.n_genes_per_strain, n_bites_i))
        genes_2 = reshape(s.host_genes[:, gene_inds_2[:]], (p.n_loci, p.n_genes_per_strain, n_bites_i))
    
        # Recombined strain IDs.
        # Reuse IDs for paired strains with the same ID. Make new IDs for the other pairs.
        # Dimensions: (n_bites_i,)
        same_strain = strain_ids_1 .== strain_ids_2
        not_same_strain = .!same_strain
        n_same = sum(same_strain)
        strain_ids_recombined =
            (same_strain .* strain_ids_1) .+
            fill_mask_with_entries(not_same_strain, next_strain_ids!(s, n_same))
    
        # Recombined genes.
        # For non-recombining infections (given by same_strain), shuffle by column.
        # For recombining infections, sample from genes_1 and genes_2 stacked on top of each other.
        genes_recombined = fill(0, (p.n_loci, p.n_genes_per_strain, n_bites_i))
        genes_recombined[:, :, same_strain] = shuffle_columns(genes_1[:, :, same_strain])
        genes_recombined[:, :, not_same_strain] = sample_each_column_without_replacement(
            cat(genes_1[:, :, not_same_strain], genes_2[:, :, not_same_strain]; dims = 3),
            p.n_genes_per_strain
        )
    
        dst_transmit_indices = findfirst_each_column(s.strain_id_liver[:,dst_hosts_i] .== 0)
        inf_host_indices = zip_cartesian(dst_transmit_indices, dst_hosts_i)
        s.strain_id_liver[inf_host_indices] = strain_ids_recombined
        s.genes_liver[:, inf_host_indices] = genes_recombined
    end
    
    verify(p, s)
    
    (n_infected_bites, n_bites, n_bites_raw)
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
        function (i)
            ind = genes_cartesian[i] # Gets multidimensional index from linear index
            host = ind[1]
            inf = ind[2]
            s.expression_index[host, inf] > 0
        end,
        mut_indices_raw
    )
    n_mut = length(mut_indices)
#     println("n_mut_raw = $(n_mut_raw), n_mut = $(n_mut)")

    # Apply mutations by creating new alleles
    # TODO: data-parallelize? (not helpful at low mutation rates.)
    for i in mut_indices
        (host, inf, exp_ind, locus) = Tuple(genes_cartesian[i])
        @assert s.expression_index[host, inf] > 0
        @assert maximum(s.n_alleles) < typemax(AlleleId)
        s.n_alleles[locus] += 1
        s.genes_active[host, inf, exp_ind, locus] = s.n_alleles[locus]

        s.strain_id_active[host, inf] = s.next_strain_id
        s.next_strain_id += 1
    end
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
        (host, inf) = Tuple(inf_cartesian[recomb_indices[i]])

        i1, i2 = samplepair(p.n_genes_per_strain)
        src_gene_1 = @view s.genes_active[host, inf, i1, :]
        src_gene_2 = @view s.genes_active[host, inf, i2, :]
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

            s.strain_id_active[host, inf] = s.next_strain_id
            s.next_strain_id += 1
        end
    end
end

function do_immigration!(t, p, s)
#     println("do_immigration!()")

    # TODO: adjust for dt
    p_bite = 1 - exp(
        -p.immigration_rate_fraction * p.biting_rate[1 + Int(floor(t)) % p.t_year]
    )
#     println("p_bite = $(p_bite)")
    n_bites = rand(Binomial(p.n_hosts, p_bite))


    # Sample source hosts
    hosts_raw = sample(1:p.n_hosts, n_bites, replace = false)

    # Filter to hosts with available infection slots
    hosts = filter(
        i -> sum(isnan.(s.t_infection_liver[i,:])) > 0,
        hosts_raw
    )

#     println("n_imm_raw = $(length(hosts_raw)), n_imm = $(length(hosts))")

    for host in hosts
        inf_ind = findfirst(isnan.(s.t_infection_liver[host, :]))

        s.strain_id_liver[host, inf_ind] = s.next_strain_id
        s.next_strain_id += 1
        s.genes_liver[host, inf_ind, :, :] = s.gene_pool[rand(1:size(s.gene_pool)[1], p.n_genes_per_strain), :]
        s.t_infection_liver[host, inf_ind] = t
    end
end
