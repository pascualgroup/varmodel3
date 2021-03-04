using SQLite: DB, Stmt
import SQLite.DBInterface.execute

struct Database
    db::DB
    summary::Stmt
    gene_strain_counts::Stmt
end

function execute(db::Database, cmd)
    execute(db.db, cmd)
end

function initialize_database(p::Params)
    if isfile(p.output_db_filename)
        error("$(p.output_db_filename) already exists; delete first")
    end
    db = DB(p.output_db_filename)
    
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
    
    Database(
        db,
        make_insert_statement(db, "summary", 11),
        make_insert_statement(db, "gene_strain_counts", 3)
    )
end

function make_insert_statement(db, table_name, n_columns)
    qmarks = join(repeat(["?"], n_columns), ",")
    Stmt(
        db,
        "INSERT INTO $(table_name) VALUES ($(qmarks))"
    )
end
