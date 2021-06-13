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
    sampled_infection_genes::Stmt
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
            n_circulating_genes INTEGER,
            n_circulating_strains INTEGER
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

    allele_columns = join(["allele_id_$(i) INTEGER" for i in 1:P.n_loci], ", ")
    execute(db, """
        CREATE TABLE sampled_infection_genes(
            infection_id INTEGER,
            expression_index INTEGER,
            $(allele_columns)
        );
    """)

    VarModelDB(
        db,
        make_insert_statement(db, "meta", 2),
        make_insert_statement(db, "summary", 13),
        make_insert_statement(db, "gene_strain_counts", 3),
        make_insert_statement(db, "sampled_hosts", 6),
        make_insert_statement(db, "sampled_infections", 6),
        make_insert_statement(db, "sampled_infection_genes", 2 + P.n_loci)
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


### GENERIC TOP-LEVEL OUTPUT FUNCTION ###

"""
Construct an prepared insert statement for a particular table.

The statement covers all columns in the table, and `n_columns` must match the
actual number of columns in the table.
"""
function write_output!(db, t, s, stats)
    if P.t_burnin !== missing && t < P.t_burnin
        return
    end

    if t % minimum(
        (P.summary_period, P.host_sampling_period, P.gene_strain_count_period,)
    ) == 0
        println("t = $(t)")

#         println("write_output!($(t))")

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
