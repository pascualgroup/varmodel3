"""
This file defines types and functions for database output.

Currently, it is shared across model variants, but in the future may need to be
broken up to be variant-specific.

The `write_output!()` function defined here calls three functions defined by the
code for individual model variants:

* `write_summary()`: writes out summary data
* `write_host_samples()`: writes out sampled hosts and infections
* `write_gene_strain_counts()`: writes out gene and strain counts

There is a sophisticated Julia system called Tables.jl that would allow the code
here to be less repetitive. With the goal of not introducing too many layers
of abstraction, this code uses SQLite more directly; a future Julia-oriented
maintainer may want to convert the code.
"""

using SQLite: DB, Stmt
import SQLite.DBInterface.execute
if P.calc_summary_statistics_instead_sqlite
    using DataFrames
    using LinearAlgebra
    using Statistics
end

"""
Database type encapsulating SQLite.DB and prepared insert statements for tables.

Prepared insert statements are used for the sake of performance; they allow
SQLite to avoid parsing and compiling the statements with every inserted row.
"""

struct VarModelDB
    db::DB
    meta::Stmt
    summary::Stmt
    gene_strain_counts::Stmt
    sampled_hosts::Stmt
    sampled_infections::Stmt
    sampled_durations::Stmt
    sampled_infection_genes::Stmt
    sampled_immunity::Stmt
    summary_statistics::Stmt
end


"""
Type encapsulating various summary statistics gathered between summary periods.
"""
@with_kw mutable struct SummaryStats
    start_datetime::DateTime
    n_events::Int = 0
    n_bites::Int = 0
    n_infected_bites::Int = 0
    n_infected_bites_with_space::Int = 0
    n_transmitting_bites::Int = 0
    n_transmissions::Int = 0
end

"""
SummaryStats constructor.
"""
function SummaryStats()
    SummaryStats(start_datetime = now())
end

"""
Reset counts and time in a SummaryStats object.

Called after summary statistics are written out in model variant-specific
`write_summary()` function.
"""
function reset!(stats::SummaryStats, start_datetime)
    stats.start_datetime = start_datetime
    stats.n_bites = 0
    stats.n_infected_bites = 0
    stats.n_infected_bites_with_space = 0
    stats.n_transmitting_bites = 0
    stats.n_transmissions = 0
end

"""
Pass commands issued to a `VarModelDB` on to the underlying `SQLite.DB`.
"""
function execute(db::VarModelDB, cmd)
    execute(db.db, cmd)
end

"""
Initialize database.

Note: if the number of columns in a table is modified, the corresponding
`make_insert_statement()` call must be updated to match the number of columns.
(This could be automated if it seems worth it.)
"""
function initialize_database()
    if isfile(P.output_db_filename)
        error("$(P.output_db_filename) already exists; delete first")
    end
    db = DB(P.output_db_filename)

    # Note: all of this could be done with the Julia Tables package, which
    # might make it less cumbersome, but also introduces conceptual overhead. 
       
    execute(db, "CREATE TABLE meta (key, value);")

    execute(db, """
        CREATE TABLE summary (
            time INTEGER,
            n_infections_liver INTEGER,
            n_infections_active INTEGER,
            n_infections INTEGER,
            n_infected_liver INTEGER,
            n_infected_active INTEGER,
            n_infected INTEGER,
            n_bites INTEGER,
            n_infected_bites INTEGER,
            n_infected_bites_with_space INTEGER,
            n_transmitting_bites INTEGER,
            n_transmissions INTEGER,
            exec_time INTEGER
        );
    """)

    execute(db, """
        CREATE TABLE gene_strain_counts (
            time INTEGER,
            n_circulating_genes_liver INTEGER,
            n_circulating_strains_liver INTEGER,
            n_circulating_genes_blood INTEGER,
            n_circulating_strains_blood INTEGER
        );
    """)

    execute(db, """
        CREATE TABLE sampled_hosts (
            time INTEGER,
            id INTEGER,
            birth_time REAL,
            death_time REAL,
            n_infections_liver INTEGER,
            n_infections_active INTEGER
        )
    """)

    execute(db, """
        CREATE TABLE sampled_infections (
            time INTEGER,
            host_id INTEGER,
            infection_id INTEGER,
            infection_time REAL,
            strain_id INTEGER,
            expression_index INTEGER
        );
    """)

    execute(db, """
        CREATE TABLE sampled_durations (
            host_id INTEGER,
            n_cleared_infections INTEGER,
            n_immune_alleles INTEGER,
            infection_time REAL,
            expression_time REAL,
            infection_duration REAL
        );
    """)

    allele_columns = join(["allele_id_$(i) INTEGER" for i in 1:P.n_loci], ", ")
    execute(db, """
        CREATE TABLE sampled_infection_genes(
            infection_id INTEGER,
            expression_index INTEGER,
            group_id INTEGER,
            $(allele_columns),
            UNIQUE(infection_id, expression_index) ON CONFLICT IGNORE
        );
    """)

    execute(db, """
        CREATE TABLE sampled_immunity (
            time INTEGER,
            host_id INTEGER,
            locus INTEGER,
            allele INTEGER,
            immunity_level INTEGER
        );
    """)

    execute(db, """
        CREATE TABLE summary_statistics (
            time INTEGER,
            prevalence REAL,
            meanMOIvar REAL,
            meanPTS REAL,
            meanPTSGroupA REAL,
            meanPTSGroupBC REAL,
            numStrains INTEGER,
            numGenes INTEGER,
            numGenesGroupA INTEGER,
            numGenesGroupBC INTEGER
        );
    """)

    VarModelDB(
        db,
        make_insert_statement(db, "meta", 2),
        make_insert_statement(db, "summary", 13),
        make_insert_statement(db, "gene_strain_counts", 5),
        make_insert_statement(db, "sampled_hosts", 6),
        make_insert_statement(db, "sampled_infections", 6),
        make_insert_statement(db, "sampled_durations", 6),
        make_insert_statement(db, "sampled_infection_genes", 2 + 1 + P.n_loci),
        make_insert_statement(db, "sampled_immunity", 5),
        make_insert_statement(db, "summary_statistics", 10)
    )
end

"""
Construct an prepared insert statement for a particular table.

The statement covers all columns in the table, and `n_columns` must match the
actual number of columns in the table.
"""
function make_insert_statement(db, table_name, n_columns)
    qmarks = join(repeat(["?"], n_columns), ",")
    Stmt(
        db,
        "INSERT INTO $(table_name) VALUES ($(qmarks))"
    )
end


### OUTPUT FUNCTIONS ###

function write_output!(db, t, s, stats)
    if P.t_burnin !== nothing && t < P.t_burnin
        return
    end

    if P.calc_summary_statistics_instead_sqlite 
        if t in P.calc_summary_statistics_times
            println(stderr, "t = $(t)")

            execute(db, "BEGIN TRANSACTION")

            write_summary_statistics(db, t, s)

            execute(db, "COMMIT")
            flush(stderr)
            flush(stdout)
        end
        return
    end

    if t % P.summary_period == 0 || t % P.gene_strain_count_period == 0 || (((P.t_host_sampling_start !== nothing && t >= P.t_host_sampling_start) || P.t_host_sampling_start === nothing) && t % P.t_year in P.host_sampling_period)
        println(stderr, "t = $(t)")

        execute(db, "BEGIN TRANSACTION")

        if P.generalized_immunity_on && t % P.summary_period == 0 && ((P.t_host_sampling_start !== nothing && t >= P.t_host_sampling_start) || P.t_host_sampling_start === nothing) && t % P.t_year in P.host_sampling_period
            hosts_GI_impact = Dict(host.id => Int[] for host in s.hosts)
            for host in s.hosts
                if length(host.active_infections) > 0
                    GI_impact_vector = rand(s.rng, Float64, length(host.active_infections)) .< exp(-host.generalized_immunity * P.generalized_immunity_detectability_param)
                    hosts_GI_impact[host.id] = Int.(GI_impact_vector)
                end
            end
        end        

        if t % P.summary_period == 0
            if P.generalized_immunity_on && ((P.t_host_sampling_start !== nothing && t >= P.t_host_sampling_start) || P.t_host_sampling_start === nothing) && t % P.t_year in P.host_sampling_period 
                write_summary(db, t, s, stats, hosts_GI_impact)
            else
                write_summary(db, t, s, stats)
            end
            write_duration!(db, t, s)
        end

        if t % P.gene_strain_count_period  == 0 
            if P.generalized_immunity_on && ((P.t_host_sampling_start !== nothing && t >= P.t_host_sampling_start) || P.t_host_sampling_start === nothing) && t % P.t_year in P.host_sampling_period 
                write_gene_strain_counts(db, t, s, hosts_GI_impact)
            else
                write_gene_strain_counts(db, t, s)
            end
        end    
        
        if ((P.t_host_sampling_start !== nothing && t >= P.t_host_sampling_start) || P.t_host_sampling_start === nothing) && t % P.t_year in P.host_sampling_period
            if P.generalized_immunity_on
                write_host_samples(db, t, s, hosts_GI_impact)
            else
                write_host_samples(db, t, s)
            end
        end        

        execute(db, "COMMIT")
        flush(stderr)
        flush(stdout)
    end
end

"""
Write output to `summary` table
"""
function write_summary(db, t, s, stats)

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
Write output to `summary` table, with generalized immunity on
"""
function write_summary(db, t, s, stats, hosts_GI_impact)

    # Compute number of infections (liver, active, and both) with a simple tally/sum.
    n_infections_liver = sum(length(host.liver_infections) for host in s.hosts)
    n_infections_active = sum(sum(hosts_GI_impact[host.id]) for host in s.hosts)
    n_infections = n_infections_liver + n_infections_active

    # Compute number of individuals with infections (liver, active, or either).
    n_infected_liver = sum(length(host.liver_infections) > 0 for host in s.hosts)
    n_infected_active = sum(any(hosts_GI_impact[host.id] >. 0) for host in s.hosts)
    n_infected = sum(
        length(host.liver_infections) > 0 || any(hosts_GI_impact[host.id] >. 0)
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
Write output to `summary_statistics` table
"""
function write_summary_statistics(db, t, s)
    sampled_hosts = sample(s.hosts, P.host_sample_size, replace = false)
    sampled_hosts_infected = filter(x -> length(x.active_infections) > 0, sampled_hosts)
    if length(sampled_hosts_infected) == 0
        return false
    end
    if P.generalized_immunity_on
        sampled_hosts_infected_GI_impact = Dict(shi.id => Int[] for shi in sampled_hosts_infected)
        sampled_hosts_infected_GI = []
        for shi in sampled_hosts_infected
            GI_impact_vector = rand(s.rng, Float64, length(shi.active_infections)) .< exp(-shi.generalized_immunity * P.generalized_immunity_detectability_param)
            sampled_hosts_infected_GI_impact[shi.id] = Int.(GI_impact_vector)
            if any(GI_impact_vector .> 0)
                push!(sampled_hosts_infected_GI, shi)
            end
        end
    else 
        sampled_hosts_infected_GI = sampled_hosts_infected
    end
    # println(length(sampled_hosts_infected))
    if length(sampled_hosts_infected_GI) == 0
        return false
    end
    if P.p_microscopy_detection == nothing
        sampled_hosts_infected_GI_detected = sampled_hosts_infected_GI
    else
        sampled_hosts_infected_GI_detected = filter(x -> rand(s.rng) < P.p_microscopy_detection, sampled_hosts_infected_GI)
    end
    # println(length(sampled_hosts_infected_detected))
    if length(sampled_hosts_infected_GI_detected) == 0
        return false
    end
    # prevalence
    preval = length(sampled_hosts_infected_GI_detected)/length(sampled_hosts)

    # number of strains and genes for blood-stage infections
    sampled_infections_detected = DataFrame(host_id = Int[], infection_id = Int[], strain_id = Int[], index = Int[], gene_id = String[], group_id = Int[])
    if P.generalized_immunity_on
        if P.undersampling_of_var
            for host in sampled_hosts_infected_GI_detected
                GI_impact_vector = sampled_hosts_infected_GI_impact[host.id]
                for infection_index in 1:length(host.active_infections)
                    infection = host.active_infections[infection_index]
                    if GI_impact_vector[infection_index] > 0
                        write_infection_df_measurement_error!(sampled_infections_detected, host, infection, s) 
                    end
                end
            end
        else 
            for host in sampled_hosts_infected_GI_detected
                GI_impact_vector = sampled_hosts_infected_GI_impact[host.id]
                for infection_index in 1:length(host.active_infections)
                    infection = host.active_infections[infection_index]
                    if GI_impact_vector[infection_index] > 0
                        write_infection_df!(sampled_infections_detected, host, infection, s)
                    end
                end
            end
        end
    else 
        if P.undersampling_of_var
            for host in sampled_hosts_infected_GI_detected
                for infection in host.active_infections
                    write_infection_df_measurement_error!(sampled_infections_detected, host, infection, s)
                end 
            end
        else 
            for host in sampled_hosts_infected_GI_detected
                for infection in host.active_infections
                    write_infection_df!(sampled_infections_detected, host, infection, s)
                end
            end
        end
    end
    sampled_infections_detected_groupA = sampled_infections_detected[isequal.(sampled_infections_detected.group_id, 1),:]
    sampled_infections_detected_groupBC = sampled_infections_detected[isequal.(sampled_infections_detected.group_id, 2),:]

    nbgene = length(unique(sampled_infections_detected[:,"gene_id"]))
    nbgeneA = length(unique(sampled_infections_detected_groupA[:,"gene_id"]))
    nbgeneBC = length(unique(sampled_infections_detected_groupBC[:,"gene_id"]))
    nbstrain = length(unique(sampled_infections_detected[:,"strain_id"]))
    
    MOI = DataFrame(HostID = Int[], MOI = Int[], Prob = Float64[])
    for host in sampled_hosts_infected_detected
        # println(host.id)
        sampled_infections_detected_host = sampled_infections_detected_groupBC[isequal.(sampled_infections_detected_groupBC.host_id, host.id),:]
        # println(unique(sampled_infections_detected_host[:,"host_id"]))
        isolateSize = length(unique(sampled_infections_detected_host[:,"gene_id"]))
        # println(isolateSize)
        moi = calcMOI(isolateSize, host, P.MOI_aggregate_approach)
        append!(MOI, moi)
    end
    if P.MOI_aggregate_approach == "pool"
        meanMOIvar = mean(MOI[:, "MOI"])
    else
        meanMOIvar = sum(MOI[:, "MOI"].*MOI[:, "Prob"])/length(unique(MOI[:, "HostID"]))
    end

    sampled_infections_detected_duplicatedGenesRemoved = sampled_infections_detected[.!nonunique(DataFrame(sampled_infections_detected[:,["host_id", "gene_id"]])),:]
    PTS = calcPTS(sampled_infections_detected_duplicatedGenesRemoved)
    sampled_infections_detected_groupA_duplicatedGenesRemoved = sampled_infections_detected_groupA[.!nonunique(DataFrame(sampled_infections_detected_groupA[:,["host_id", "gene_id"]])),:]
    PTSGroupA = calcPTS(sampled_infections_detected_groupA_duplicatedGenesRemoved)
    sampled_infections_detected_groupBC_duplicatedGenesRemoved = sampled_infections_detected_groupBC[.!nonunique(DataFrame(sampled_infections_detected_groupBC[:,["host_id", "gene_id"]])),:]
    PTSGroupBC = calcPTS(sampled_infections_detected_groupBC_duplicatedGenesRemoved)
    meanPTS = mean(offdiag(PTS))
    meanPTSGroupA = mean(offdiag(PTSGroupA))
    meanPTSGroupBC = mean(offdiag(PTSGroupBC))

    execute(db.summary_statistics, (
        t,
        preval,
        meanMOIvar,
        meanPTS,
        meanPTSGroupA,
        meanPTSGroupBC,
        nbstrain,
        nbgene,
        nbgeneA,
        nbgeneBC
    ))
end

function write_infection_df_measurement_error!(df, host, infection, s)
    unique_genes_A = Set()
    unique_genes_BC = Set()
    for i in 1:size(infection.genes)[2]
        gene_group_id = s.association_genes_to_var_groups[Gene(infection.genes[:,i])]
        if gene_group_id == 1 
            push!(unique_genes_A, infection.genes[:,i])
        else
            push!(unique_genes_BC, infection.genes[:,i])
        end
    end
    num_A_detected = rand(s.rng, P.measurement_error_A)
    while num_A_detected > length(unique_genes_A)
        num_A_detected = rand(s.rng, P.measurement_error_A)
    end 
    unique_genes_A_detected = sample(collect(unique_genes_A), num_A_detected, replace = false)
    gene_group_id = 1
    for i in 1:length(unique_genes_A_detected)
        # println(unique_genes_A_detected[i])
        gene_id_rename = join(unique_genes_A_detected[i], "-")
        # println(gene_id_rename)
        index = i
        push!(df, [host.id infection.id infection.strain_id index gene_id_rename gene_group_id])
    end
    num_BC_detected = rand(s.rng, P.measurement_error_BC)
    while num_BC_detected > length(unique_genes_BC)
        num_BC_detected = rand(s.rng, P.measurement_error_BC)
    end 
    unique_genes_BC_detected = sample(collect(unique_genes_BC), num_BC_detected, replace = false)
    gene_group_id = 2
    for j in 1:length(unique_genes_BC_detected)
        gene_id_rename = join(unique_genes_BC_detected[j], "-")
        index = length(unique_genes_A_detected) + j
        push!(df, [host.id infection.id infection.strain_id index gene_id_rename gene_group_id])
    end
end

function write_infection_df!(df, host, infection, s)
    unique_genes = Set()
    for i in 1:size(infection.genes)[2]
        push!(unique_genes, infection.genes[:,i])
    end
    for i in 1:length(unique_genes)
        gene_id = Gene(unique_genes[i])
        # @assert haskey(s.association_genes_to_var_groups, gene_id)
        gene_group_id = s.association_genes_to_var_groups[gene_id]
        gene_id_rename = join(unique_genes[i], "-") 
        index = i
        push!(df, [host.id infection.id infection.strain_id index gene_id_rename gene_group_id])
    end
end

function calcMOI(isolateSize, host, MOI_aggregate_approach)
    denominator = 0.0
    numerators = []
    for i in 1:P.maxMOI 
        prob1_dict = P.p_isolateSize_given_MOI[i]
        if string(isolateSize) in keys(prob1_dict)
            prob1 = prob1_dict[string(isolateSize)]
        else
            prob1 = 0.0
        end
        prob2 = P.MOI_prior[i]
        denominator_temp = prob1 * prob2
        push!(numerators, denominator_temp)
        denominator += denominator_temp
    end
    probMOIGivenIsoSize = numerators./denominator
    if MOI_aggregate_approach == "pool"
        moi = DataFrame(HostID = host.id, MOI = argmax(probMOIGivenIsoSize), Prob = maximum(probMOIGivenIsoSize))
    else
        moi = DataFrame(HostID = host.id, MOI = 1:P.maxMOI, Prob = probMOIGivenIsoSize)
    end
    moi
end

function calcPTS(infections)
    df1 = unstack(infections, :host_id, :gene_id, :index) # https://dataframes.juliadata.org/stable/man/reshaping_and_pivoting/
    df2 = select(df1, Not(:host_id))
    df3 = coalesce.(df2, 0)
    df4 = df3.>0
    df5 = df4.*1
    df6 = Matrix(df5)
    overlap = df6 * transpose(df6)
    isolateSize = sum(df6, dims = 2)
    pts = overlap./isolateSize[:,1]
    pts
end

function offdiag(A::Matrix)
    @assert size(A)[1] == size(A)[2]
    D = size(A)[1]
    v = zeros(D*(D-1))
    for i in 1:D
        for j in 1:(i-1)
            v[(i-1)*(D-1)+j] = A[j,i]
        end
        for j in (i+1):D
            v[(i-1)*(D-1)+j-1] = A[j,i]
        end
    end
    v
end

"""
Write output for periodically sampled hosts.
"""
function write_host_samples(db, t, s)

    # Sample `host_sample_size` hosts randomly (without replacement).
    hosts = sample(s.hosts, P.host_sample_size, replace = false)

    # For each host, write out birth/death time and each infection.
    for host in hosts
        write_immunity(db, t, host)

        execute(
        db.sampled_hosts,
        (
            t, Int64(host.id), host.t_birth, missing, # host.t_death, # Cannot know t_death
            length(host.liver_infections), length(host.active_infections)
        )
        )

        for infection in host.liver_infections
            write_infection(db, t, host, infection, s)
        end

        for infection in host.active_infections
            write_infection(db, t, host, infection, s)
        end

    end
end


"""
Write output for periodically sampled hosts, with generalized immunity on.
"""
function write_host_samples(db, t, s, hosts_GI_impact)

    # Sample `host_sample_size` hosts randomly (without replacement).
    hosts = sample(s.hosts, P.host_sample_size, replace = false)

    # For each host, write out birth/death time and each infection.
    for host in hosts
        host_GI_impact = hosts_GI_impact[host.id]
        if any(host_GI_impact .> 0)
            write_immunity(db, t, host)

            execute(
            db.sampled_hosts,
            (
                t, Int64(host.id), host.t_birth, missing, # host.t_death, # Cannot know t_death
                length(host.liver_infections), sum(host_GI_impact)
            )
            )

            for infection in host.liver_infections
                write_infection(db, t, host, infection, s)
            end

            for infection_index in 1:length(host.active_infections)
                if host_GI_impact[infection_index] > 0
                    infection = host.active_infections[infection_index]
                    write_infection(db, t, host, infection, s)
                end
            end
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
function write_infection(db, t, host, infection, s)
    execute(
        db.sampled_infections,
        (
            t, Int64(host.id), Int64(infection.id), infection.t_infection, Int64(infection.strain_id),
            infection.expression_index == 0 ? missing : Int64(infection.expression_index)
        )
    )
    for i in 1:P.n_genes_per_strain
        gene_id = Gene(infection.genes[:,i])
        # @assert haskey(s.association_genes_to_var_groups, gene_id)
        gene_group_id = s.association_genes_to_var_groups[gene_id]
        execute(db.sampled_infection_genes, vcat([infection.id, i], gene_group_id, infection.genes[:,i]))
    end
end

"""
Write the infection durations to the `sampled_durations` table.
"""
function write_duration!(db, t, s)
    for dur in s.durations
        execute(
            db.sampled_durations,
            (
                Int64(dur.host_id), Int64(dur.n_cleared_infections),
                Int64(dur.n_immune_alleles),
                dur.t_infection, dur.t_expression, dur.duration
            )
        )
    end
    empty!(s.durations)
end

"""
Write gene and strain counts to the `gene_strain_counts` table.

Counts are not maintained dynamically during the simulation; this function
simply scans all host infections and assembles sets of all genes and all
strains.
"""
function write_gene_strain_counts(db, t, s)
    genesLiver::Set{Gene} = Set()
    strainsLiver::BitSet = BitSet()
    genesBlood::Set{Gene} = Set()
    strainsBlood::BitSet = BitSet()

    for host in s.hosts
        count_genes_and_strains!(genesLiver, strainsLiver, host.liver_infections)
        count_genes_and_strains!(genesBlood, strainsBlood, host.active_infections)
    end

    execute(db.gene_strain_counts, (t, length(genesLiver), length(strainsLiver), length(genesBlood), length(strainsBlood)))
end

"""
Write gene and strain counts to the `gene_strain_counts` table.

Counts are not maintained dynamically during the simulation; this function
simply scans all host infections and assembles sets of all genes and all
strains.

Generalized immunity is on.
"""
function write_gene_strain_counts(db, t, s, hosts_GI_impact)
    genesLiver::Set{Gene} = Set()
    strainsLiver::BitSet = BitSet()
    genesBlood::Set{Gene} = Set()
    strainsBlood::BitSet = BitSet()

    for host in s.hosts
        count_genes_and_strains!(genesLiver, strainsLiver, host.liver_infections)
        host_GI_impact = hosts_GI_impact[host.id]
        for infection_index in 1:length(host.active_infections)
            if host_GI_impact[infection_index] > 0
                count_genes_and_strains!(genesBlood, strainsBlood, host.active_infections[infection_index])
            end
        end
    end

    execute(db.gene_strain_counts, (t, length(genesLiver), length(strainsLiver), length(genesBlood), length(strainsBlood)))
end


"""
Count genes and strains for a particular list of infections in a particular host.

This function simply adds genes and strains from each infection to corresponding
sets.
"""
function count_genes_and_strains!(genes, strains, infections)
    for infection in infections
        for i in 1:P.n_genes_per_strain
            push!(genes, infection.genes[:,i])
        end
        push!(strains, infection.strain_id)
    end
end

"""
Write immunity.
"""
function write_immunity(db, t, host)
    for locus in 1:P.n_loci
        d = host.immunity.vd[locus]
        for (key, value) in d
            execute(db.sampled_immunity, (t, Int64(host.id), Int64(locus), Int64(key), Int64(value)))
        end
    end
end
