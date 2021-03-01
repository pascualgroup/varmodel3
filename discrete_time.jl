function run_discrete_time(p::Params)
    s = State(p)
    db = initialize_database(p)
    
    execute(db, "BEGIN TRANSACTION")
    write_summary(0, p, s, db)
    execute(db, "COMMIT")
    
    verify(p, s)
    
    execute(db, "BEGIN TRANSACTION")
    
    # TODO: use dt
    for t in 1:p.t_end
#         println("stepping to t = $(t)")
        
        do_rebirth!(t, p, s)
        do_immunity_loss!(t, p, s) 
        do_transition!(t, p, s)
        do_biting!(t, p, s)
        do_mutation!(t, p, s)
        do_recombination!(t, p, s)
        do_immigration!(t, p, s)
        
        if t % 360 == 0
            println("")
            println("t = $(t)")
            println("n_infections = $(sum(s.expression_index .> 0))")
            println("n_immunity = $(sum(s.immunity .> 0))")
            
            write_summary(t, p, s, db)
            
            execute(db, "COMMIT")
            
            verify(p, s)
            
            execute(db, "BEGIN TRANSACTION")
        end
    end
end

function write_summary(t, p, s::State, db)
    n_infections = sum(s.expression_index .!== 0)
    n_infected = sum(sum(s.expression_index .!== 0, dims = 2) .!== 0)
    n_infected_bites = 0
    n_total_bites = 0
    (n_circulating_genes, n_circulating_strains) = count_circulating_genes_and_strains(p, s)
    println("n_genes = $(n_circulating_genes), n_strains = $(n_circulating_strains)")
    exec_time = 0.0
    
    execute(db.summary, [
        t,
        n_infections,
        n_infected,
        n_infected_bites,
        n_total_bites,
        n_circulating_strains,
        n_circulating_genes,
        exec_time
    ])
end

function count_circulating_genes_and_strains(p::Params, s::State)
    println("Starting count...")
    genes::Set{SVector{p.n_loci, AlleleId}} = Set()
    strains::BitSet = BitSet()
    
    for j in 1:p.max_infection_count
        for i in 1:p.n_hosts
            strain_id = s.infection_strain_id[i,j]
            if strain_id != 0 && !(strain_id in strains)
                for k in 1:p.n_genes_per_strain
                    push!(genes, @view(s.infection_genes[i, j, k, :]))
                end
                push!(strains, strain_id)
            end
        end
    end
    println("Count finished.")
    
    (length(genes), length(strains))
end

function do_rebirth!(t, p, s)
#     println("do_rebirth!()")
    
    dead_hosts = findall(s.t_death .< t)
#     println("n dead: $(length(dead_hosts))")
    
    s.t_birth[dead_hosts] = s.t_death[dead_hosts]
    s.t_death[dead_hosts] = s.t_birth[dead_hosts] + [draw_host_lifetime(p) for i in 1:length(dead_hosts)]
    s.t_infection[dead_hosts, :, :] .= NaN32
    s.infection_strain_id[dead_hosts, :] .= 0
    s.infection_genes[dead_hosts, :, :, :] .= 0
    s.expression_index[dead_hosts, :] .= 0
    s.immunity[dead_hosts, :, :] .= 0
end

function do_immunity_loss!(t, p, s)
    # TODO: adjust for dt
    p_loss = 1 - exp(-p.immunity_loss_rate)
    
    # Lose immunity one locus at a time
    for locus in 1:p.n_loci
        n_immunities = p.n_hosts * s.n_alleles[locus]
        n_loss = rand(Binomial(n_immunities, p_loss))
#         println("n_loss($(locus)) = $(n_loss)")
        
        # Uniformly randomly sample indices in one-dimensional array
        indices = sample(1:n_immunities, n_loss, replace = false)
        
        # Create view of immunities for this locus as a one-dimensional array
        immunity_view = reshape(
            (@view(s.immunity[:, 1:s.n_alleles[locus], locus])),
            n_immunities
        )
        
        # Decrement immunity (leaving zeros at zero)
        immunity_view[indices] = max(zeros(ImmunityCount, n_loss), immunity_view[indices] .- ImmunityCount(1))
    end
end

function findfirst_stable(x)
    y = findfirst(x)
    y == nothing ? 0 : y
end

function do_transition!(t, p, s)
#     println("do_transition!()")
    # TODO: adjust for dt
    p_transition = 1 - exp(-p.transition_rate)
    
    n_infections = p.n_hosts * p.max_infection_count
    
    n_trans_raw = rand(Binomial(n_infections, p_transition))
#     println("n_trans_raw: $(n_trans_raw)")
    if n_trans_raw == 0
        return
    end
    
    # Uniformly randomly sample indices in one-dimensional array
    # Includes nonexistent infections
    indices_raw = sample(1:n_infections, n_trans_raw, replace = false)
    
    # Get one-dimensional views on infection
    expression_index_view = reshape(s.expression_index, n_infections)
    t_infection_view = reshape(s.t_infection, n_infections)
    infection_genes_view = reshape(s.infection_genes, (n_infections, p.n_genes_per_strain, p.n_loci))
    
    # Filter indices down to infections that exist and are past the liver stage
    indices = indices_raw[
        (expression_index_view[indices_raw] .> 0) .& (t_infection_view[indices_raw] .+ p.t_liver_stage .< t)
    ]
    n_trans = length(indices)
#     println("n_trans = $(n_trans)")
    if n_trans == 0
        return
    end
    
    # Compute host and infection indices from sampled linear indices
    host_indices = ((indices .- 1) .% p.n_hosts) .+ 1
    inf_indices = ((indices .- 1) .÷ p.n_hosts) .+ 1
#     host_inf_indices = [CartesianIndex(x) for x in zip(host_indices, inf_indices)]
    
    # Advance expression for each transition
    # TODO: make CUDA-friendly
    max_immunity_count = ImmunityCount(p.max_immunity_count)
    for i in 1:n_trans
        host_index = host_indices[i]
        inf_index = inf_indices[i]
        
        # Advance expression until non-immune gene is reached
        while true
            exp_index = s.expression_index[host_index, inf_index]
            @assert exp_index > 0
            
            # Gain extra immunity to current gene
            for locus in 1:p.n_loci
                allele_id = s.infection_genes[host_index, inf_index, exp_index, locus]
                s.immunity[host_index, allele_id, locus] = min(
                    max_immunity_count,
                    s.immunity[host_index, allele_id, locus] + 1
                )
            end
            
            if exp_index == p.n_genes_per_strain
                # Clear infection if we're at the end
                s.t_infection[host_index, inf_index] = NaN32
                s.infection_strain_id[host_index, inf_index] = 0
                s.infection_genes[host_index, inf_index, :, :] .= 0
                s.expression_index[host_index, inf_index] = 0
                break
            else
                # Advance expression if we're not yet at the end
                s.expression_index[host_index, inf_index] += 1
                
                # If we're not immune to this gene, stop advancing expression
                if !is_immune(p, s, host_index, inf_index, exp_index + 1)
                    break
                end
            end
        end
    end
end

function is_immune(p, s, host_index, inf_index, exp_index)
    for locus in 1:p.n_loci
        allele_id = s.infection_genes[host_index, inf_index, exp_index, locus]
        if s.immunity[host_index, allele_id, locus] == 0
            return false
        end
    end
    return true
end

function do_biting!(t, p, s)
    # TODO: adjust for dt
    p_bite = 1 - exp(
        -(p.biting_rate_mean * p.daily_biting_rate_multiplier[
            1 + Int(floor(t)) % p.t_year
        ])
    )
#     println("p_bite = $(p_bite)")
    n_bites = rand(Binomial(p.n_hosts, p_bite))
    
    # Sample source hosts and destination hosts
    src_hosts_raw = sample(1:p.n_hosts, n_bites, replace = false)
    dst_hosts_raw = sample(1:p.n_hosts, n_bites, replace = false)
    if length(src_hosts_raw) == 0
        return
    end
    
    # Identify indices with active source infections and available space in destination
    src_dst_valid = let
        t_infection = s.t_infection[src_hosts_raw, :]
        src_valid = fill(false, n_bites)
        for i in 1:p.max_infection_count
            src_valid .|= (t_infection[:, i] .+ p.t_liver_stage .< t)
        end
        
        expression_index = s.expression_index[dst_hosts_raw, :]
        dst_valid = fill(false, n_bites)
        for i in 1:p.max_infection_count
            dst_valid .|= expression_index[:, i] .== 0
        end
        
        src_valid .& dst_valid
    end
    src_hosts = src_hosts_raw[src_dst_valid]
    if length(src_hosts) == 0
        return
    end
    dst_hosts = dst_hosts_raw[src_dst_valid]
    
    # Compute number of active source infections and their locations in array
    src_inf_count, src_inf_indices = let
        t_infection = s.t_infection[src_hosts, :]
        count = fill(0, length(src_hosts))
        is_active = fill(false, length(src_hosts))
        
        next_index = fill(1, length(src_hosts))
        indices = fill(0, (length(src_hosts), p.max_infection_count))
        for i in 1:p.max_infection_count
            is_active .= t_infection[:, i] .+ p.t_liver_stage .< t
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
    
    # Compute number of available infection slots in destination and their locations in array
    dst_inf_count, dst_inf_indices = let
        expression_index = s.expression_index[dst_hosts, :]
        count = fill(0, length(dst_hosts))
        is_available = fill(false, length(dst_hosts))
        
        next_index = fill(1, length(dst_hosts))
        indices = fill(0, (length(src_hosts), p.max_infection_count))
        for i in 1:p.max_infection_count
            is_available .= expression_index[:, i] .== 0
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
            
            strain1 = reshape(@view(s.infection_genes[src_host, i1, :, :]), p.n_genes_per_strain, p.n_loci)
            strain2 = reshape(@view(s.infection_genes[src_host, i2, :, :]), p.n_genes_per_strain, p.n_loci)
            if s.infection_strain_id[src_host, i1] == s.infection_strain_id[src_host, i2]
                # Just copy shuffled genes if the strains are identical
                s.infection_strain_id[dst_host, dst_inf_index] = s.infection_strain_id[src_host, i1]
                s.infection_genes[dst_host, dst_inf_index, :, :] = strain1[randperm(p.n_genes_per_strain), :]
            else
                # Recombine genes if the strains are not identical
                all_genes = vcat(strain1, strain2)
                s.infection_strain_id[dst_host, dst_inf_index] = s.next_strain_id
                s.next_strain_id += 1
                s.infection_genes[dst_host, dst_inf_index, :, :] = all_genes[sample(1:(2 * p.n_genes_per_strain), p.n_genes_per_strain, replace = true), :]
            end
            s.expression_index[dst_host, dst_inf_index] = 1
            s.t_infection[dst_host, dst_inf_index] = t
            
#             println("strain1 = $(strain1)")
#             println("strain2 = $(strain2)")
#             println("dst_strain = $(s.infection_genes[dst_host, dst_inf_index, :, :])")
        end
    end
end

function do_mutation!(t, p, s)
    # TODO: adjust for dt
    p_mut = 1 - exp(-p.mutation_rate)
    
    # Draw mutations randomly using linear indexing
    n_possible_mut = length(s.infection_genes)
    n_mut_raw = rand(Binomial(n_possible_mut, p_mut))
    mut_indices_raw = sample(1:n_possible_mut, n_mut_raw, replace = false)
    
    # Apply mutations for infections that actually exist
    genes_cartesian = CartesianIndices(s.infection_genes)
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
        s.infection_genes[host, inf, exp_ind, locus] = s.n_alleles[locus]
        
        s.infection_strain_id[host, inf] = s.next_strain_id
        s.next_strain_id += 1
    end
    
    resize_immunity_if_necessary!(p, s)
end

function do_recombination!(t, p, s)
    p_recomb = 1 - exp(
        -p.ectopic_recombination_rate * (p.n_genes_per_strain) * (p.n_genes_per_strain - 1) / 2.0
    )
    
    n_inf_raw = p.n_hosts * p.max_infection_count
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
        src_gene_1 = @view s.infection_genes[host, inf, i1, :]
        src_gene_2 = @view s.infection_genes[host, inf, i2, :]
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
            
            s.infection_strain_id[host, inf] = s.next_strain_id
            s.next_strain_id += 1
        end
    end
end

function do_immigration!(t, p, s)
    # TODO: adjust for dt
    p_bite = 1 - exp(
        -p.immigration_rate_fraction * (p.biting_rate_mean * p.daily_biting_rate_multiplier[
            1 + Int(floor(t)) % p.t_year
        ])
    )
#     println("p_bite = $(p_bite)")
    n_bites = rand(Binomial(p.n_hosts, p_bite))
    
    
    # Sample source hosts
    hosts_raw = sample(1:p.n_hosts, n_bites, replace = false)
    
    # Filter to hosts with available infection slots
    hosts = filter(
        i -> sum(s.expression_index[i,:] .== 0) > 0,
        hosts_raw
    )
    
#     println("n_imm_raw = $(length(hosts_raw)), n_imm = $(length(hosts))")
    
    for host in hosts
        inf_ind = findfirst(s.expression_index[host, :] .== 0)
        
        s.infection_strain_id[host, inf_ind] = s.next_strain_id
        s.next_strain_id += 1
        s.infection_genes[host, inf_ind, :, :] = s.gene_pool[rand(1:size(s.gene_pool)[1], p.n_genes_per_strain), :]
        s.t_infection[host, inf_ind] = t
        s.expression_index[host, inf_ind] = 1
    end
end

function resize_immunity_if_necessary!(p, s)
    immunity_size = size(s.immunity)[2]
    if maximum(s.n_alleles) > immunity_size
        println("increasing immunity capacity by 25%...")
        new_immunity_size = immunity_size * 5 ÷ 4
        new_immunity = fill(ImmunityCount(0), p.n_hosts, new_immunity_size, p.n_loci)
        new_immunity[:,1:immunity_size,:] = s.immunity
        s.immunity = new_immunity
    end
end
