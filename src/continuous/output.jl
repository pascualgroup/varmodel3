"""
Output functions for exact continuous-time model.
"""

"""
Write output to `summary` table
"""
function write_summary(db, t, s, stats)
#     println("write_summary($(t))")

    # Compute number of infections (liver, active, and both) with a simple tally/sum.
    n_infections_liver = sum(length(host.liver_infections) for host in s.hosts)
    n_infections_active = sum(length(host.active_infections) for host in s.hosts)
    n_infections = n_infections_liver + n_infections_active

    # Compute number of individuals with infections (liver, active, or either).
    n_infected_liver = sum(length(host.liver_infections) > 0 for host in s.hosts)
    n_infected_active = sum(length(host.active_infections) > 0 for host in s.hosts)
    n_infected = sum(
        length(host.liver_infections) > 0 || length(host.active_infections) > 0
        for host in s.hosts
    )

    # Compute elapsed time in seconds.
    next_datetime = now()
    exec_time = Dates.value(next_datetime - stats.start_datetime) / 1000.0

    # Write to summary table.
    execute(db.summary, (
        t,
        n_infections_liver,
        n_infections_active,
        n_infections,
        n_infected_liver,
        n_infected_active,
        n_infected,
        stats.n_bites,
        stats.n_infected_bites,
        stats.n_infected_bites_with_space,
        stats.n_transmitting_bites,
        stats.n_transmissions,
        exec_time
    ))

    # Reset counters and elapsed time.
    reset!(stats, next_datetime)
end

"""
Write output for periodically sampled hosts.
"""
function write_host_samples(db, t, s)
#     println("write_host_samples($(t))")

    # Sample `host_sample_size` hosts randomly (without replacement).
    hosts = sample(s.hosts, P.host_sample_size, replace = false)

    # For each host, write out birth/death time and each infection.
    for host in hosts
        execute(
            db.sampled_hosts,
            (
                t, Int64(host.id), host.t_birth, host.t_death,
                length(host.liver_infections), length(host.active_infections)
            )
        )

        for infection in host.liver_infections
            write_infection(db, t, host, infection)
        end

        for infection in host.active_infections
            write_infection(db, t, host, infection)
        end
    end
end

"""
Write output for a single infection from a sampled host.

Infections are written to the `sampled_infections` table, and the genes for the
infections are written to `sampled_infection_genes`.

`expression_index` indicates which gene is currently expressed. Liver-stage
infections are indicated using a SQLite `NULL` (Julia `missing`) for
`expression_index`.

For each infection, the `sampled_infection_genes` table contains one row for
each expression index, with the final columns containing the allele IDs for each
locus.
"""
function write_infection(db, t, host, infection)
#     println("write_infection($(t), $(host.id), $(infection.id))")
    execute(
        db.sampled_infections,
        (
            t, Int64(host.id), Int64(infection.id), infection.t_infection, Int64(infection.strain_id),
            infection.expression_index == 0 ? missing : Int64(infection.expression_index)
        )
    )
    for i in 1:P.n_genes_per_strain
        execute(db.sampled_infection_genes, vcat([infection.id, i], infection.genes[:,i]))
    end
end

"""
Write gene and strain counts to the `gene_strain_counts` table.

Counts are not maintained dynamically during the simulation; this function
simply scans all host infections and assembles sets of all genes and all
strains.
"""
function write_gene_strain_counts(db, t, s)
#     println("write_gene_strain_counts($(t))")
    genes::Set{Gene} = Set()
    strains::BitSet = BitSet()

    for host in s.hosts
        count_genes_and_strains!(genes, strains, host.liver_infections)
        count_genes_and_strains!(genes, strains, host.active_infections)
    end

    execute(db.gene_strain_counts, (t, length(genes), length(strains)))
end

"""
Count genes and strains for a particular list of infections in a particular host.

This function simply adds genes and strains from each infection to corresponding
sets.
"""
function count_genes_and_strains!(genes, strains, infections)
#     println("count_circulating_genes_and_strains()")
    for infection in infections
        for i in 1:P.n_genes_per_strain
            push!(genes, infection.genes[:,i])
        end
        push!(strains, infection.strain_id)
    end
end
