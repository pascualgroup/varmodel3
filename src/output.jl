using SQLite: DB, Stmt
import SQLite.DBInterface.execute

struct Database
    db::DB
    meta::Stmt
    summary::Stmt
    gene_strain_counts::Stmt
    sampled_hosts::Stmt
    sampled_infections::Stmt
    sampled_infection_genes::Stmt
end

@with_kw mutable struct SummaryStats
    start_datetime::DateTime
    n_events::Int = 0
    n_bites::Int = 0
    n_infected_bites::Int = 0
    n_infected_bites_with_space::Int = 0
    n_transmitting_bites::Int = 0
    n_transmissions::Int = 0
end

function SummaryStats()
    SummaryStats(start_datetime = now())
end

function reset!(stats::SummaryStats, start_datetime)
    stats.start_datetime = start_datetime
    stats.n_bites = 0
    stats.n_infected_bites = 0
    stats.n_infected_bites_with_space = 0
    stats.n_transmitting_bites = 0
    stats.n_transmissions = 0
end

function execute(db::Database, cmd)
    execute(db.db, cmd)
end

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
            n_circulating_genes INTEGER,
            n_circulating_strains INTEGER
        );
    """)
    
    execute(db, """
        CREATE TABLE sampled_hosts (
            time INTEGER,
            id INTEGER,
            birth_time REAL,
            death_time REAL
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
    
    allele_columns = join(["allele_id_$(i) INTEGER" for i in 1:P.n_loci], ", ")
    execute(db, """
        CREATE TABLE sampled_infection_genes(
            infection_id INTEGER,
            expression_index INTEGER,
            $(allele_columns)
        );
    """)
    
    Database(
        db,
        make_insert_statement(db, "meta", 2),
        make_insert_statement(db, "summary", 13),
        make_insert_statement(db, "gene_strain_counts", 3),
        make_insert_statement(db, "sampled_hosts", 4),
        make_insert_statement(db, "sampled_infections", 6),
        make_insert_statement(db, "sampled_infection_genes", 2 + P.n_loci)
    )
end

function make_insert_statement(db, table_name, n_columns)
    qmarks = join(repeat(["?"], n_columns), ",")
    Stmt(
        db,
        "INSERT INTO $(table_name) VALUES ($(qmarks))"
    )
end


### GENERIC TOP-LEVEL OUTPUT FUNCTION ###

function write_output(db, t, s, stats)
    if t % minimum(
        (P.summary_period, P.host_sampling_period, P.gene_strain_count_period,)
    ) == 0
        println("t = $(t)")
        
#         println("write_output($(t))")
        
        execute(db, "BEGIN TRANSACTION")
        
        if t % P.summary_period == 0
            write_summary(db, t, s, stats)
        end
        
        if t % P.host_sampling_period == 0
            write_host_samples(db, t, s)
        end
        
        if t % P.gene_strain_count_period == 0
            write_gene_strain_counts(db, t, s)
        end
        
        execute(db, "COMMIT")
    end
end
