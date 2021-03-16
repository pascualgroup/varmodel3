using SQLite: DB, Stmt
import SQLite.DBInterface.execute
import SQLite.sqlreturn


sqlreturn(context, val::Int8) = sqlreturn(context, Int64(val))
sqlreturn(context, val::UInt8) = sqlreturn(context, Int64(val))
sqlreturn(context, val::Int16) = sqlreturn(context, Int64(val))
sqlreturn(context, val::UInt16) = sqlreturn(context, Int64(val))
sqlreturn(context, val::UInt32) = sqlreturn(context, Int64(val))
sqlreturn(context, val::UInt64) = sqlreturn(context, Int64(val))
# sqlreturn(context, val::T) where {T <: Integer} = sqlreturn(context, Int64(val))

struct Database
    db::DB
    summary::Stmt
    gene_strain_counts::Stmt
    sampled_hosts::Stmt
    sampled_infections::Stmt
    sampled_infection_genes::Stmt
end

function execute(db::Database, cmd)
    execute(db.db, cmd)
end

function initialize_database(p::Params)
    if isfile(p.output_db_filename)
        error("$(p.output_db_filename) already exists; delete first")
    end
    db = DB(p.output_db_filename)
    
    # Note: all of this could be done with the Julia Tables package, which
    # might make it less cumbersome, but also introduces conceptual overhead.
    
    execute(db, "CREATE TABLE meta (key, value);")
    
    execute(db, """
        CREATE TABLE summary (
            time REAL,
            n_infections_liver INTEGER,
            n_infections_active INTEGER,
            n_infections INTEGER,
            n_infected_liver INTEGER,
            n_infected_active INTEGER,
            n_infected INTEGER,
            n_infected_bites INTEGER,
            n_infected_bites_with_space INTEGER,
            n_total_bites INTEGER,
            exec_time INTEGER
        );
    """)
    
    execute(db, """
        CREATE TABLE gene_strain_counts (
            time REAL,
            n_circulating_genes INTEGER,
            n_circulating_strains INTEGER
        );
    """)
    
    execute(db, """
        CREATE TABLE sampled_hosts (
            time REAL,
            id INTEGER,
            birth_time REAL,
            death_time REAL
        )
    """)
    
    execute(db, """
        CREATE TABLE sampled_infections (
            time REAL,
            host_id INTEGER,
            infection_id INTEGER,
            infection_time REAL,
            strain_id INTEGER,
            expression_index INTEGER
        );
    """)
    
    allele_columns = join(["allele_id_$(i) INTEGER" for i in 1:p.n_loci], ", ")
    execute(db, """
        CREATE TABLE sampled_infection_genes(
            infection_id INTEGER,
            expression_index INTEGER,
            $(allele_columns)
        );
    """)
    
    Database(
        db,
        make_insert_statement(db, "summary", 11),
        make_insert_statement(db, "gene_strain_counts", 3),
        make_insert_statement(db, "sampled_hosts", 4),
        make_insert_statement(db, "sampled_infections", 6),
        make_insert_statement(db, "sampled_infection_genes", 2 + p.n_loci)
    )
end

function make_insert_statement(db, table_name, n_columns)
    qmarks = join(repeat(["?"], n_columns), ",")
    Stmt(
        db,
        "INSERT INTO $(table_name) VALUES ($(qmarks))"
    )
end
