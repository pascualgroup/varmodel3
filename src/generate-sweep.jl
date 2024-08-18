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
import DataStructures.OrderedDict
using DataFrames

# Get relevant paths and cd to the script path.
# NB: use actual relative locations of varmodel3 root relative to your script.
include("../../preamble.jl")
SCRIPT_PATH = abspath(dirname(PROGRAM_FILE))
ROOT_PATH = abspath(joinpath(SCRIPT_PATH, "..", ".."))
ROOT_RUN_SCRIPT = joinpath(ROOT_PATH, "run.jl")
ROOT_RUNMANY_SCRIPT = joinpath(ROOT_PATH, "runmany.jl")
cd(SCRIPT_PATH)

# Number of replicates for each parameter combination.
const N_REPLICATES = parse(Int, ARGS[3])
PARAMSTXT_DIR = ARGS[1]
PARAMSTXT_FILE = ARGS[2]
# NosRun = eval(parse(ARGS[4]))    
NosRun = let expr = Meta.parse(ARGS[4])
    @assert expr.head == :vect
    Int.(expr.args)
end
# println(NosRun)

# Number of SLURM jobs to generate.
const N_JOBS_MAX = size(NosRun)[1]
const N_CORES_PER_JOB_MAX = N_REPLICATES # Half a node, easier to get scheduled than a whole one

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
    # sqlite does not have a separate boolean storage class, store as integer 0/1 instead, more efficient than true/false string.
    execute(db, "CREATE TABLE param_combos (combo_id INTEGER, No INTEGER, daily_biting_rate_multiplier_file TEXT, irs_start_year INTEGER,
                                            irs_duration INTEGER, biting_rate_factor REAL, t_end_years INTEGER,
                                            t_host_sampling_start_year INTEGER, n_genes_initial INTEGER,
                                            n_genes_per_strain INTEGER, n_alleles_per_locus_initial INTEGER,
                                            var_groups_do_not_share_alleles INTEGER, var_groups_ratio_A REAL, var_groups_ratio_BC REAL,
                                            var_groups_fix_ratio INTEGER, var_groups_functionality_A REAL, var_groups_functionality_BC REAL,
                                            ectopic_recombination_rate_A REAL, ectopic_recombination_rate_BC REAL, ectopic_recombination_generates_new_alleles INTEGER,
                                            p_ectopic_recombination_generates_new_allele REAL, var_groups_high_functionality_express_earlier INTGER,
                                            biting_rate_mean REAL, immigration_rate_fraction REAL, switching_rate_A REAL, switching_rate_BC REAL)") 
    execute(db, "CREATE TABLE runs (run_id INTEGER, combo_id INTEGER, replicate INTEGER, rng_seed INTEGER, run_dir TEXT, params TEXT)")
    execute(db, "CREATE TABLE jobs (job_id INTEGER, job_dir TEXT)")
    execute(db, "CREATE TABLE job_runs (job_id INTEGER, run_id INTEGER)")

    generate_runs(db)
    generate_jobs(db)
end

function generate_runs(db)
    # System random device used to generate seeds.
    seed_rng = RandomDevice()
    
    paramstxtFull_dir = PARAMSTXT_DIR
    paramstxtFull_file = PARAMSTXT_FILE
    paramstxtFull_path = joinpath(paramstxtFull_dir, paramstxtFull_file)
    paramstxtFull_df, paramstxtFull_header = readdlm(paramstxtFull_path, '\t', header = true)
    paramstxtFull = DataFrame(paramstxtFull_df, vec(paramstxtFull_header))
    paramstxt = paramstxtFull[ [x in NosRun for x in paramstxtFull[:,"No"]],:]
    
    # Base parameter set, copied/modified for each combination/replicate.
    base_params = init_base_params()
    validate(base_params)
    execute(db, "INSERT INTO meta VALUES (?, ?)", ("base_params", pretty_json(base_params)))

    # Loop through parameter combinations and replicates, generating a run directory
    # `runs/c<combo_id>/r<replicate>` for each one.
    combo_id = 1
    run_id = 1
    for row in 1:size(paramstxt, 1)
        No = paramstxt[row, "No"]
        daily_biting_rate_multiplier_file = paramstxt[row, "daily_biting_rate_multiplier_file"]
        irs_start_year = paramstxt[row, "irs_start_year"]
        irs_duration = paramstxt[row, "irs_duration"]
        biting_rate_factor = paramstxt[row, "biting_rate_factor"]
        t_end_years = paramstxt[row, "t_end_years"]
        t_host_sampling_start_year = paramstxt[row, "t_host_sampling_start_year"]
        host_sample_size = paramstxt[row, "host_sample_size"]
        n_genes_initial = paramstxt[row, "n_genes_initial"]
        n_genes_per_strain = paramstxt[row, "n_genes_per_strain"]
        n_alleles_per_locus_initial = paramstxt[row, "n_alleles_per_locus_initial"]
        var_groups_do_not_share_alleles = paramstxt[row, "var_groups_do_not_share_alleles"] == "True"
        var_groups_ratio_A = paramstxt[row, "var_groups_ratio_A"] 
        var_groups_ratio_BC = paramstxt[row, "var_groups_ratio_BC"] 
        var_groups_ratio = [var_groups_ratio_A, var_groups_ratio_BC]
        var_groups_fix_ratio = paramstxt[row, "var_groups_fix_ratio"]=="True"
        var_groups_functionality_A = paramstxt[row, "var_groups_functionality_A"]
        var_groups_functionality_BC = paramstxt[row, "var_groups_functionality_BC"]
        var_groups_functionality = [var_groups_functionality_A, var_groups_functionality_BC]
        ectopic_recombination_rate_A = paramstxt[row, "ectopic_recombination_rate_A"]
        ectopic_recombination_rate_BC = paramstxt[row, "ectopic_recombination_rate_BC"]
        ectopic_recombination_rate = [ectopic_recombination_rate_A, ectopic_recombination_rate_BC]
        ectopic_recombination_generates_new_alleles = paramstxt[row, "ectopic_recombination_generates_new_alleles"] == "True"
        if ectopic_recombination_generates_new_alleles
            p_ectopic_recombination_generates_new_allele = paramstxt[row, "p_ectopic_recombination_generates_new_allele"]
        else
            p_ectopic_recombination_generates_new_allele = nothing
        end
        var_groups_high_functionality_express_earlier = paramstxt[row, "var_groups_high_functionality_express_earlier"] == "True"
        biting_rate_mean = paramstxt[row, "biting_rate_mean"]
        immigration_rate_fraction = paramstxt[row, "immigration_rate_fraction"]
        t_year = paramstxt[row, "t_year"]
        switching_rate_A = paramstxt[row, "switching_rate_A"]
        switching_rate_BC = paramstxt[row, "switching_rate_BC"]
        switching_rate = [switching_rate_A, switching_rate_BC]
        gene_strain_count_period = paramstxt[row, "gene_strain_count_period"]
        n_initial_infections = paramstxt[row, "n_initial_infections"]
        t_burnin_years = paramstxt[row, "t_burnin_years"]
        migrants_match_local_prevalence = paramstxt[row, "migrants_match_local_prevalence"] == "True"
        upper_bound_recomputation_period = paramstxt[row, "upper_bound_recomputation_period"]
        immunity_loss_rate = paramstxt[row, "immunity_loss_rate"]
        mutation_rate = paramstxt[row, "mutation_rate"]
        liver_erlang_shape = paramstxt[row, "liver_erlang_shape"]

        t_end = t_end_years * t_year
        t_host_sampling_start = t_host_sampling_start_year * t_year

        if isfile(daily_biting_rate_multiplier_file) 
            daily_biting_rate_multiplier = readdlm(daily_biting_rate_multiplier_file, Float64)[:,1]
        else
            daily_biting_rate_multiplier = repeat([1.0], t_year)
        end

        biting_rate = biting_rate_mean * daily_biting_rate_multiplier

        biting_rate_multiplier_by_year = repeat([1.0], t_end_years)
        if irs_start_year !== nothing
            index = (irs_start_year + 1):(irs_start_year + irs_duration)
            biting_rate_multiplier_by_year[index] .= biting_rate_factor
        end

        host_sampling_period = [180, 300]

        t_burnin = t_burnin_years * t_year
        
        
        calc_summary_statistics_instead_sqlite = paramstxt[row, "calc_summary_statistics_instead_sqlite"] == "True"
        calc_summary_statistics_times = readdlm(paramstxt[row, "calc_summary_statistics_times_file"], Int64)[:,1]
        MOI_aggregate_approach = paramstxt[row, "MOI_aggregate_approach"]
        maxMOI = paramstxt[row, "maxMOI"]
        MOI_prior = readdlm(paramstxt[row, "MOI_prior_file"], Float64)[:,1]
        p_isolateSize_given_MOI = []
        for i in 1:maxMOI
            p_isolateSize_given_MOI_df_file = string(paramstxt[row, "p_isolateSize_given_MOI_file"]) * string(i) * ".txt"
            p_isolateSize_given_MOI_df = readdlm(p_isolateSize_given_MOI_df_file)
            isolateSizes = Array{Int64}(p_isolateSize_given_MOI_df[:,1])
            Probs = Array{Float64}(p_isolateSize_given_MOI_df[:,2])
            # p_isolateSize_given_MOI_dict = Dict(k => Float64(0) for k in isolateSizes)
            p_isolateSize_given_MOI_dict = Dict{String, Float64}()
            for j in 1:length(isolateSizes)
                isolateSize = isolateSizes[j]
                p_isolateSize_given_MOI_dict[string(isolateSize)] = Probs[j]
            end
            push!(p_isolateSize_given_MOI, p_isolateSize_given_MOI_dict)
        end
        p_microscopy_detection = paramstxt[row, "p_microscopy_detection"]
        undersampling_of_var = paramstxt[row, "undersampling_of_var"] == "True"
        measurement_error_A = readdlm(paramstxt[row, "measurement_error_A_file"], Int64)[:,1]
        measurement_error_BC = readdlm(paramstxt[row, "measurement_error_BC_file"], Int64)[:,1]

        println("Processing c$(combo_id): No = $(No)")
        execute(db, "INSERT INTO param_combos VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (combo_id, No, daily_biting_rate_multiplier_file, irs_start_year, irs_duration, biting_rate_factor, t_end_years, t_host_sampling_start_year, n_genes_initial, n_genes_per_strain, n_alleles_per_locus_initial, var_groups_do_not_share_alleles, var_groups_ratio_A, var_groups_ratio_BC, var_groups_fix_ratio, var_groups_functionality_A, var_groups_functionality_BC, ectopic_recombination_rate_A, ectopic_recombination_rate_BC, ectopic_recombination_generates_new_alleles, p_ectopic_recombination_generates_new_allele, var_groups_high_functionality_express_earlier, biting_rate_mean, immigration_rate_fraction, switching_rate_A, switching_rate_BC))

        for replicate in 1:N_REPLICATES
            output_db_filename = "sim_" * string(No) * "_r" * string(replicate) * "_sd.sqlite" 
            rng_seed = rand(seed_rng, 1:typemax(Int64))
            params = add_params(base_params, (
                rng_seed = rng_seed,
                t_year = t_year,
                t_end = t_end,
                t_host_sampling_start = t_host_sampling_start,
                host_sample_size = host_sample_size,
                t_burnin = t_burnin,
                upper_bound_recomputation_period = upper_bound_recomputation_period,

                biting_rate = biting_rate, 
                irs_start_year = irs_start_year,
                irs_duration = irs_duration,
                biting_rate_multiplier_by_year = biting_rate_multiplier_by_year, 
                
                n_genes_initial = n_genes_initial,
                n_genes_per_strain = n_genes_per_strain,
                n_alleles_per_locus_initial = n_alleles_per_locus_initial,

                var_groups_do_not_share_alleles = var_groups_do_not_share_alleles,
                var_groups_ratio = var_groups_ratio,
                var_groups_fix_ratio = var_groups_fix_ratio,
                var_groups_functionality = var_groups_functionality,
                ectopic_recombination_rate = ectopic_recombination_rate,
                ectopic_recombination_generates_new_alleles = ectopic_recombination_generates_new_alleles,
                p_ectopic_recombination_generates_new_allele = p_ectopic_recombination_generates_new_allele,
                var_groups_high_functionality_express_earlier = var_groups_high_functionality_express_earlier,
                switching_rate = switching_rate,

                immigration_rate_fraction = immigration_rate_fraction,

                gene_strain_count_period = gene_strain_count_period,
                host_sampling_period = host_sampling_period,

                output_db_filename = output_db_filename,

                n_initial_infections = n_initial_infections,
                migrants_match_local_prevalence = migrants_match_local_prevalence,
                immunity_loss_rate = immunity_loss_rate,
                mutation_rate = mutation_rate,
                liver_erlang_shape = liver_erlang_shape,
                
                calc_summary_statistics_instead_sqlite = calc_summary_statistics_instead_sqlite,
                calc_summary_statistics_times = calc_summary_statistics_times,
                MOI_aggregate_approach = MOI_aggregate_approach,
                maxMOI = maxMOI,
                MOI_prior = MOI_prior,
                p_isolateSize_given_MOI = p_isolateSize_given_MOI,
                p_microscopy_detection = p_microscopy_detection,
                undersampling_of_var = undersampling_of_var,
                measurement_error_A = measurement_error_A,
                measurement_error_BC = measurement_error_BC
            ))

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

            #SBATCH --account=pi-jozik
            #SBATCH --partition=broadwl

            #SBATCH --job-name=var-$(job_id)

            #SBATCH --tasks=1
            #SBATCH --cpus-per-task=$(n_cores)
            #SBATCH --mem-per-cpu=6000m
            #SBATCH --time=36:00:00

            #SBATCH --chdir=$(joinpath(SCRIPT_PATH, job_dir))
            #SBATCH --output=output.txt

            # Uncomment this to use the Midway-provided Julia:
            module load julia/1.7.2

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
    d = OrderedDict(fn => getfield(params, fn) for fn in fieldnames(typeof(params)))
    io = IOBuffer()
    JSON.print(io, d, 2)
    String(take!(io))
end

function init_base_params()
    biting_rate_mean = 0.00005
    daily_biting_rate_multiplier = readdlm("../mosquito_population.txt", Float64)[:,1]
    biting_rate = biting_rate_mean * daily_biting_rate_multiplier
    t_year = 360 
    t_burnin_year = 2
    t_burnin = t_burnin_year * t_year
    t_host_sampling_start_year = 197
    t_host_sampling_start = t_host_sampling_start_year * t_year
    t_end_years = 200
    t_end = t_end_years * t_year
    biting_rate_factor = 1.0
    biting_rate_multiplier_by_year =  repeat([biting_rate_factor], t_end_years)
    n_genes_initial = 20000

    add_params(Params(), ( 
        t_year = t_year,
        t_end = t_end,
        t_burnin = t_burnin,
        t_host_sampling_start = t_host_sampling_start,

        upper_bound_recomputation_period = 30,

        output_db_filename = "output.sqlite",

        summary_period = 30,
        gene_strain_count_period = 30,

        host_sampling_period = [180,300],
        host_sample_size = 1500,

        verification_period = t_end/10,

        sample_infection_duration_every = 100,

        rng_seed = nothing,

        whole_gene_immune = false,

        n_hosts = 10000,
        n_initial_infections = 20,

        n_genes_initial = n_genes_initial,
        n_genes_per_strain = 60,

        n_loci = 2,

        n_alleles_per_locus_initial = n_genes_initial/10,

        transmissibility = 1.0,
        coinfection_reduces_transmission = true,

        ectopic_recombination_rate = [0.000742, 0.000742],
        p_ectopic_recombination_is_conversion = 0.0,

        ectopic_recombination_generates_new_alleles = false,
        p_ectopic_recombination_generates_new_allele = nothing,

        rho_recombination_tolerance = 0.8,
        mean_n_mutations_per_epitope = 5.0,

        immunity_loss_rate = 1/1080,

        mutation_rate = 1.42e-8,

        t_liver_stage = 14.0,
        liver_erlang_shape = 49,

        switching_rate = [1.0/6.0, 1.0/6.0],

        mean_host_lifetime = 23.38837487739662 * t_year,

        background_clearance_rate = 0.0,

        immigration_rate_fraction = 0.0,

        n_infections_liver_max = 20,
        n_infections_active_max = 20,
        
        biting_rate = biting_rate,

        irs_start_year = 201,
        irs_duration = 0,
        biting_rate_multiplier_by_year = biting_rate_multiplier_by_year,

        migrants_match_local_prevalence = true,
        migration_rate_update_period = 30,

        # parameters for var groups implementation
        var_groups_functionality = [1.0,1.0],
        var_groups_ratio = [1.0,0.0],
        var_groups_fix_ratio = false,
        var_groups_do_not_share_alleles = false,
        var_groups_high_functionality_express_earlier = false,
        gene_group_id_association_recomputation_period = 30,
        
        # parameters for fitting
        calc_summary_statistics_instead_sqlite = false,
        calc_summary_statistics_times = nothing,
        MOI_aggregate_approach = nothing,
        maxMOI = nothing,
        MOI_prior = nothing,
        p_isolateSize_given_MOI = nothing,
        p_microscopy_detection = nothing,
        undersampling_of_var = false,
        measurement_error_A = nothing,
        measurement_error_BC = nothing,
    ))
end

main()