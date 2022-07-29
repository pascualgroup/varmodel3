#!/usr/bin/env julia

"""
The purpose of this file is to illustrate how to do a parameter sweep in Julia.

This script loops through parameter combinations, and replicates with different
random seeds, and generates files necessary to perform runs on a local machine
or on a SLURM cluster.

To use this script for an experiment, you should copy this directory to a new
location, modify the parameter sweeps, and modify the relative paths to
`preamble.jl` and `ROOT_PATH` below to be correct, and the run

`./generate-runs.jl`

When the experiment is complete, you can collect all output files into a single
SQLite database via

`./gather-output.jl`

For each run, it creates a directory, `runs/c<combo_id>/r<replicate>`, and adds
entries to a SQLite database of run information, to make it easy to identify
runs and collate output.

It also divides runs into jobs suitable for execution on a single cluster node
or local machine. The runs are specified as lines in the job's `runs.txt`
file, and the job is specified in a `job.sbatch` file, which can be run directly
as a shell script or submitted to a SLURM cluster.
Each job uses the script `varmodel3/runmany.jl` to run a single-node, multi-core
queue of runs, with one run running on each core at any time.

This script also generates a script `submit_jobs.sh`, which submits every job to
SLURM at once.

Runs are divided into at most `N_JOBS_MAX` jobs that make use of at most
`N_CORES_PER_JOB_MAX` for the cluster node's local queue.
This allows you to work within limits set by your cluster administrator.
If you have no limits, you should set `N_JOBS_MAX` to a very large number,
and set `N_CORES_PER_JOB_MAX = 1`, so that the cluster can dynamically
balance runs across cluster nodes as the experiment runs.

To modify configuration settings for SLURM jobs, edit the template string in
the `generate_jobs()` function.
"""

println("(Annoying Julia compilation delay...)")

using Random
using SQLite
import SQLite.DBInterface.execute
using DelimitedFiles

# Get relevant paths and cd to the script path.
# NB: use actual relative locations of varmodel3 root relative to your script.
include("../../preamble.jl")
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_PATH = abspath(joinpath(SCRIPT_PATH, "..", ".."))
ROOT_RUN_SCRIPT = joinpath(ROOT_PATH, "run.jl")
ROOT_RUNMANY_SCRIPT = joinpath(ROOT_PATH, "runmany.jl")
cd(SCRIPT_PATH)

# Number of replicates for each parameter combination.
const N_REPLICATES = 2

# Number of SLURM jobs to generate.
const N_JOBS_MAX = 100
const N_CORES_PER_JOB_MAX = 14 # Half a node, easier to get scheduled than a whole one

function main()
    # Root run directory.
    if ispath("runs")
        error("`runs` already exists; please move or delete.")
    end
    mkdir("runs")

    # Root job directory.
    if ispath("jobs")
        error("`jobs` already exists; please move or delete.")
    end
    mkdir("jobs")

    # Database of experiment information.
    if ispath("sweep_db.sqlite")
        error("`sweep_db.sqlite` already exists; please move or delete")
    end
    db = SQLite.DB(joinpath("sweep_db.sqlite"))
    execute(db, "CREATE TABLE meta (key, value)")
    execute(db, "CREATE TABLE param_combos (combo_id INTEGER, mutation_rate REAL, transmissibility REAL)")
    execute(db, "CREATE TABLE runs (run_id INTEGER, combo_id INTEGER, replicate INTEGER, rng_seed INTEGER, run_dir TEXT, params TEXT)")
    execute(db, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(db, "CREATE TABLE job_runs (job_id INTEGER, run_id INTEGER)")

    generate_runs(db)
    generate_jobs(db)
end

function generate_runs(db)
    # System random device used to generate seeds.
    seed_rng = RandomDevice()

    # Base parameter set, copied/modified for each combination/replicate.
    base_params = init_base_params()
    validate(base_params)
    execute(db, "INSERT INTO meta VALUES (?, ?)", ("base_params", pretty_json(base_params)))

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    combo_id = 1
    run_id = 1
    for mutation_rate in (0.5e-8, 1.0e-8, 1.5e-8, 2.0e-8)
        for transmissibility in (0.25, 0.5, 0.75)
            println("Processing c$(combo_id): mutation_rate = $(mutation_rate), transmissibility = $(transmissibility)")

            execute(db, "INSERT INTO param_combos VALUES (?, ?, ?)", (combo_id, mutation_rate, transmissibility))

            for replicate in 1:N_REPLICATES
                rng_seed = rand(seed_rng, 1:typemax(Int64))
                params = Params(
                    base_params;
                    rng_seed = rng_seed,
                    mutation_rate = mutation_rate,
                    transmissibility = transmissibility
                )

                run_dir = joinpath("runs", "c$(combo_id)", "r$(replicate)")
                @assert !ispath(run_dir)
                mkpath(run_dir)

                # Generate parameters file.
                params_json = pretty_json(params)
                open(joinpath(run_dir, "parameters.json"), "w") do f
                    println(f, params_json)
                end

                # Generate shell script to perform a single run.
                run_script = joinpath(run_dir, "run.sh")
                open(run_script, "w") do f
                    print(f, """
                    #!/bin/sh

                    cd `dirname \$0`
                    julia --check-bounds=no -O3 $(ROOT_RUN_SCRIPT) parameters.json &> output.txt
                    """)
                end
                run(`chmod +x $(run_script)`) # Make run script executable

                # Save all run info (including redundant stuff for reference) into DB.
                execute(db, "INSERT INTO runs VALUES (?, ?, ?, ?, ?, ?)", (run_id, combo_id, replicate, rng_seed, run_dir, params_json))

                run_id += 1
            end
            combo_id += 1
        end
    end
end

function generate_jobs(db)
    println("Assigning runs to jobs...")

    # Assign runs to jobs (round-robin).
    job_id = 1
    for (run_id, run_dir) in execute(db, "SELECT run_id, run_dir FROM runs ORDER BY replicate, combo_id")
        execute(db, "INSERT INTO job_runs VALUES (?,?)", (job_id, run_id))

        # Mod-increment job ID.
        job_id = (job_id % N_JOBS_MAX) + 1
    end

    # Create job directories containing job scripts and script to submit all jobs.
    submit_file = open("submit_jobs.sh", "w")
    println(submit_file, """
    #!/bin/sh

    cd `dirname \$0`
    """)
    for (job_id,) in execute(db, "SELECT DISTINCT job_id FROM job_runs ORDER BY job_id")
        job_dir = joinpath("jobs", "$(job_id)")
        mkpath(job_dir)

        # Get all run directories for this job.
        run_dirs = [run_dir for (run_dir,) in execute(db,
            """
            SELECT run_dir FROM job_runs, runs
            WHERE job_runs.job_id = ?
            AND runs.run_id = job_runs.run_id
            """,
            (job_id,)
        )]
        n_cores = min(length(run_dirs), N_CORES_PER_JOB_MAX)

        # Write out list of runs.
        open(joinpath(job_dir, "runs.txt"), "w") do f
            for run_dir in run_dirs
                run_script = joinpath(SCRIPT_PATH, run_dir, "run.sh")
                println(f, run_script)
            end
        end

        # Create job sbatch file.
        job_sbatch = joinpath(job_dir, "job.sbatch")
        open(job_sbatch, "w") do f
            print(f, """
            #!/bin/sh

            #SBATCH --account=pi-pascualmm
            #SBATCH --partition=broadwl

            #SBATCH --job-name=var-$(job_id)

            #SBATCH --tasks=1
            #SBATCH --cpus-per-task=$(n_cores)
            #SBATCH --mem-per-cpu=2000m
            #SBATCH --time=4:00:00

            #SBATCH --chdir=$(joinpath(SCRIPT_PATH, job_dir))
            #SBATCH --output=output.txt

            # Uncomment this to use the Midway-provided Julia:
            # module load julia

            julia $(ROOT_RUNMANY_SCRIPT) $(n_cores) runs.txt
            """)
        end
        run(`chmod +x $(job_sbatch)`) # Make run script executable (for local testing)

        execute(db, "INSERT INTO jobs VALUES (?,?)", (job_id, job_dir,))
        println(submit_file, "sbatch $(job_sbatch)")
    end
    close(submit_file)
    run(`chmod +x submit_jobs.sh`) # Make submit script executable
end

function pretty_json(params)
    d = Dict(fn => getfield(params, fn) for fn in fieldnames(typeof(params)))
    io = IOBuffer()
    JSON.print(io, d, 2)
    String(take!(io))
end

function init_base_params()
    t_year = 360
    daily_biting_rate_multiplier = readdlm("../mosquito_population.txt", Float64)[:,1]
    #snp_ld_matrix = readdlm("../pairwise_ld_coefficient_24snps.txt", Float64)

    Params(
        upper_bound_recomputation_period = 30,

        output_db_filename = "output.sqlite",

        summary_period = 30,
        gene_strain_count_period = 360,

        host_sampling_period = 30,
        host_sample_size = 100,

        verification_period = 360,

        sample_duration = 1000,

        rng_seed = nothing,
        use_immunity_by_allele = false,

#         use_immunity_by_allele = true,
#         whole_gene_immune = false,

        t_year = t_year,
        t_end = (111) * t_year,

        t_burnin = 61 * t_year,

        n_hosts = 10000,
        n_initial_infections = 20,

        n_genes_initial = 9600,
        n_genes_per_strain = 60,

        n_loci = 2,

        n_alleles_per_locus_initial = 960,

        transmissibility = 0.5,
        coinfection_reduces_transmission = true,

        ectopic_recombination_rate = 1.8e-7,
        p_ectopic_recombination_is_conversion = 0.0,

        ectopic_recombination_generates_new_alleles = false,

#         ectopic_recombination_generates_new_alleles = true,
#         p_ectopic_recombination_generates_new_allele = 0.5,

        rho_recombination_tolerance = 0.8,
        mean_n_mutations_per_epitope = 5.0,

        immunity_level_max = 100,
        immunity_loss_rate = 0.001,

        mutation_rate = 1.42e-8,

        t_liver_stage = 14.0,

        switching_rate = 1.0/6.0,

        mean_host_lifetime = 30 * t_year,
        max_host_lifetime = 80 * t_year,

        immigration_rate_fraction = 0.0026,

        n_infections_liver_max = 10,
        n_infections_active_max = 10,

        biting_rate = 0.0005 * daily_biting_rate_multiplier,

        migrants_match_local_prevalence = true,
        migration_rate_update_period = 30,

        n_snps_per_strain = 24,

        distinct_initial_snp_allele_frequencies = false,
#         distinct_initial_snp_allele_frequencies = true,
#         initial_snp_allele_frequency = [0.1, 0.9],

        snp_linkage_disequilibrium = false,
#         snp_linkage_disequilibrium = true,
#         snp_pairwise_ld = snp_ld_matrix,
    )
end

main()
