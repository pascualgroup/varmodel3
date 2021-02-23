using SQLite: DB, Stmt
import SQLite.DBInterface.execute

struct Database
    db::DB
    summary::Stmt
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
            n_infections INTEGER,
            n_infected INTEGER,
            n_infected_bites INTEGER,
            n_total_bites INTEGER,
            n_circulating_strains INTEGER,
            n_circulating_genes INTEGER,
            exec_time REAL
        );
    """)
    
    Database(
        db,
        make_insert_statement(db, "summary", 8)
    )
end

function make_insert_statement(db, table_name, n_columns)
    qmarks = join(repeat(["?"], n_columns), ",")
    Stmt(
        db,
        "INSERT INTO $(table_name) VALUES ($(qmarks))"
    )
end
