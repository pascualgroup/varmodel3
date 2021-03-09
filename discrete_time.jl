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
    s.t_infection_liver[dead_hosts, :, :] .= NaN32
    s.t_infection_active[dead_hosts, :, :] .= NaN32
    s.strain_id_liver[dead_hosts, :] .= 0
    s.strain_id_active[dead_hosts, :] .= 0
    s.genes_liver[dead_hosts, :, :, :] .= 0
    s.genes_active[dead_hosts, :, :, :] .= 0
    s.expression_index[dead_hosts, :] .= 0
    # TODO: reset immunity
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
    
    # TODO: abstract out index-matching algorithm for reuse?
    # TODO: clean up that horrible subindexing?
    
    # Identify hosts with infections ready to become active
    ready_mask_raw = s.t_infection_liver .+ p.t_liver_stage .<= t
    host_indices = findall(reshape(sum(ready_mask_raw, dims = 2), :) .> 0)
#     println("size(host_indices) = $(size(host_indices))")
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
#     println("n_hosts = $(length(host_indices))")
#     println("n_activations = $(length(src_indices))")
#     println("n_activations = $(length(dst_indices))")
    
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
#     println("n_active = $(sum(s.expression_index))")
    
    strain_id_active[dst_indices] = strain_id_liver[src_indices]
    genes_active[dst_indices, :, :] = genes_liver[src_indices, :, :]
    
    # Deactivate all ready infections (including those that didn't get activated)
    t_infection_liver[src_indices] .= NaN32
    strain_id_liver[src_indices] .= 0
    genes_liver[src_indices, :, :] .= 0
    
#     expression_index = @view(s.expression_index[host_indices, :])
#     t_infection_liver = @view(s.t_infection_liver[host_indices, :])
#     t_infection_active = @view(s.t_infection_active[host_indices, :])
#     strain_id_liver = @view s.strain_id_liver[host_indices, :]
#     strain_id_active = @view s.strain_id_active[host_indices, :]
#     genes_liver = @view s.genes_liver[host_indices, :, :, :]
#     println("size(s.genes_liver) = $(size(s.genes_liver))")
#     println("size(genes_liver) = $(size(genes_liver))")
#     genes_active = @view s.genes_active[host_indices, :, :, :]
#     for i in 1:length(src_indices)
#         println("i = $(i)")
#         host, inf1 = Tuple(src_indices[i])
#         host2, inf2 = Tuple(dst_indices[i])
#         
#         println("host = $(host)")
#         println("inf1 = $(inf1)")
#         println("inf2 = $(inf2)")
#         
#         @assert host == host2
#         
#         @assert !isnan(t_infection_liver[host, inf1])
#         @assert strain_id_liver[host, inf1] != 0
#         @assert all(genes_liver[host, inf1, :, :] .> 0)
#         
#         @assert isnan(t_infection_active[host, inf2])
#         @assert strain_id_active[host, inf2] == 0
#         @assert all(genes_active[host, inf2, :, :] .== 0)
#         @assert expression_index[host, inf2] == 0
#         
#         expression_index[host, inf2] = 1
#         t_infection_active[host, inf2] = t_infection_liver[host, inf1]
#         strain_id_active[host, inf2] = strain_id_liver[host, inf1]
#         genes_active[host, inf2, :, :] = genes_liver[host, inf1, :, :]
#         
#         t_infection_liver[host, inf1] = NaN32
#         strain_id_liver[host, inf1] = 0
#         genes_liver[host, inf1, :, :] .= 0
#         println("genes_liver = $(genes_liver[host, inf1, :, :])")
#     end
    
#     println("host_indices = $(host_indices)")
#     println("src_indices = $(src_indices)")
#     
#     println("t_infection_liver = $(s.t_infection_liver[host_indices,:][src_indices])")
#     println("genes_liver = $(s.genes_liver[host_indices,:,:,:][src_indices, :, :])")
    
#     verify(p, s)
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
    infection_genes_view = reshape(s.genes_active, (n_infections, p.n_genes_per_strain, p.n_loci))
    
    # Filter indices down to active infections
    indices = indices_raw[expression_index_view[indices_raw] .> 0]
    n_trans = length(indices)
#     println("n_trans = $(n_trans)")
    if n_trans == 0
        return
    end
    
    # Compute host and infection indices from sampled linear indices
    host_indices = ((indices .- 1) .% p.n_hosts) .+ 1
    inf_indices = ((indices .- 1) .÷ p.n_hosts) .+ 1
    
    # Advance expression for each chosen infection
    # TODO: make CUDA-friendly
    for i in 1:n_trans
        host_index = host_indices[i]
        inf_index = inf_indices[i]
        
        # Advance expression until non-immune gene is reached
        while true
            exp_index = s.expression_index[host_index, inf_index]
            @assert exp_index > 0
            
            # TODO: gain extra immunity to current gene
            
            if exp_index == p.n_genes_per_strain
                # Clear infection if we're at the end
                s.t_infection_active[host_index, inf_index] = NaN32
                s.strain_id_active[host_index, inf_index] = 0
                s.genes_active[host_index, inf_index, :, :] .= 0
                s.expression_index[host_index, inf_index] = 0
                break
            else
                # Advance expression if we're not yet at the end
                s.expression_index[host_index, inf_index] += 1
                
                # TODO: If we're not immune to this gene, stop advancing expression
                # in the meantime: always stop advancing.
                break
            end
        end
    end
end

# TODO: rewrite this whole function to be cleaner and simpler
function do_biting!(t, p, s)
#     println("do_biting!()")
    
    # TODO: adjust for dt
    p_bite = 1 - exp(-p.biting_rate[1 + Int(floor(t)) % p.t_year])
#     println("p_bite = $(p_bite)")
    n_bites = rand(Binomial(p.n_hosts, p_bite))
    
    # Sample source hosts and destination hosts
    src_hosts_raw = sample(1:p.n_hosts, n_bites, replace = false)
    dst_hosts_raw = sample(1:p.n_hosts, n_bites, replace = false)
    if length(src_hosts_raw) == 0
        return (0, 0, 0)
    end
    
    # Identify indices with active source infections and available space in destination
    src_valid = reshape(any(s.expression_index[src_hosts_raw,:] .> 0, dims = 2), n_bites)
    dst_valid = reshape(any(isnan.(s.t_infection_liver[dst_hosts_raw,:]), dims = 2), n_bites)
    
    n_infected_bites = sum(src_valid)
    
    src_dst_valid = src_valid .& dst_valid
    
    src_hosts = src_hosts_raw[src_dst_valid]
    if length(src_hosts) == 0
        return (0, 0, 0)
    end
    dst_hosts = dst_hosts_raw[src_dst_valid]
    
    n_infected_bites_with_space = length(dst_hosts)
    
    # Compute number of active source infections and their locations in array
    src_inf_count, src_inf_indices = let
        expression_index = s.expression_index[src_hosts, :]
        count = fill(0, length(src_hosts))
        is_active = falses(length(src_hosts))
        
        next_index = fill(1, length(src_hosts))
        indices = fill(0, (length(src_hosts), p.n_infections_active_max))
        for i in 1:p.n_infections_active_max
            is_active .= expression_index[:,i] .> 0
            count .+= is_active
            indices[[CartesianIndex(j, next_index[j]) for j in findall(is_active)]] .= i
            next_index[is_active] .+= 1
        end
        
        count, indices
    end
    
    # Draw which source infections will be used to form transmitted strains
    src_inf_trans_count, src_inf_trans_indices = let
        p_infect = if p.coinfection_reduces_transmission
            fill(p.transmissibility, length(src_hosts))
        else
            p.transmissibility ./ src_inf_count
        end
        
        count = [rand(Binomial(src_inf_count[i], p_infect[i])) for i in 1:length(src_hosts)]
        indices = [
            [src_inf_indices[i, j] for j in shuffle(1:count[i])]
            for i in 1:length(src_hosts)
        ]
        
        (count, indices)
    end
    
    # Compute number of available infection slots in destination livers and their locations in array
    dst_inf_count, dst_inf_indices = let
        t_infection = s.t_infection_liver[dst_hosts, :]
        count = fill(0, length(dst_hosts))
        is_available = falses(length(dst_hosts))
        
        next_index = fill(1, length(dst_hosts))
        indices = fill(0, (length(src_hosts), p.n_infections_liver_max))
        for i in 1:p.n_infections_liver_max
            is_available .= isnan.(t_infection[:, i])
            count .+= is_available
            indices[[CartesianIndex(j, next_index[j]) for j in findall(is_available)]] .= i
            next_index[is_available] .+= 1
        end
        
        count, indices
    end
    
    # Compute how many infections will be transmitted
    inf_trans_count = min.(src_inf_trans_count, dst_inf_count)
    inf_trans_count_total = sum(inf_trans_count)
    
    # Construct flattened arrays whose first dimension identifies hosts
    host_subindices = let
        subindices = fill(0, inf_trans_count_total)
        k = 1
        for i in 1:length(src_hosts)
            for j in 1:inf_trans_count[i]
                subindices[k] = i
                k += 1
            end
        end
        @assert k == inf_trans_count_total + 1
        subindices
    end
    inf_src_hosts = src_hosts[host_subindices]
    inf_dst_hosts = dst_hosts[host_subindices]
    
    # Transmit:
    # TODO: make this parallel-by-host? (if slow path)
    for k in 1:inf_trans_count_total
        host_subindex = host_subindices[k]
        src_host = src_hosts[host_subindex]
        dst_host = dst_hosts[host_subindex]
        for j in 1:inf_trans_count[host_subindex]
            # Sample two infection indices from transmitting infections in source host
            (i1, i2) = rand(
                src_inf_trans_indices[host_subindex][1:src_inf_trans_count[host_subindex]], 2
            )
            
            dst_inf_index = dst_inf_indices[host_subindex, j]
            @assert dst_inf_count[host_subindex] >= j
            
            strain1 = reshape(@view(s.genes_active[src_host, i1, :, :]), p.n_genes_per_strain, p.n_loci)
            strain2 = reshape(@view(s.genes_active[src_host, i2, :, :]), p.n_genes_per_strain, p.n_loci)
            if s.strain_id_active[src_host, i1] == s.strain_id_active[src_host, i2]
                # Just copy shuffled genes if the strains are identical
                s.strain_id_liver[dst_host, dst_inf_index] = s.strain_id_active[src_host, i1]
                s.genes_liver[dst_host, dst_inf_index, :, :] = strain1[randperm(p.n_genes_per_strain), :]
            else
                # Recombine genes if the strains are not identical
                all_genes = vcat(strain1, strain2)
                s.strain_id_liver[dst_host, dst_inf_index] = s.next_strain_id
                s.next_strain_id += 1
                s.genes_liver[dst_host, dst_inf_index, :, :] = all_genes[sample(1:(2 * p.n_genes_per_strain), p.n_genes_per_strain, replace = true), :]
            end
            s.t_infection_liver[dst_host, dst_inf_index] = t
            
#             println("strain1 = $(strain1)")
#             println("strain2 = $(strain2)")
#             println("dst_strain = $(s.infection_genes[dst_host, dst_inf_index, :, :])")
        end
    end
#     println("HERE")
    (n_infected_bites, n_infected_bites_with_space, n_bites)
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