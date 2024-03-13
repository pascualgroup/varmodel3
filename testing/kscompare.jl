"""
    This file contains functions to use Kolmogorov-Smirnov tests to compare output distributions.

    The main reason to do this is to verify that output distributions are the same
    (er, not-not-the-same) between different versions of code.

    All the functions assume that output from many runs has already been gathered into
    two SQLite databases (one for each version of the code being compared).

    None of the code here is specific to varmodel3; it can be used to compare outputs between
    any SQLite databases.

    Comparisons can be made using identical queries for the two databases or different queries
    for each database, for situations where the output format is different.
"""

println("(Annoying Julia compilation delay...)")

using SQLite
import SQLite.DBInterface.execute
import HypothesisTests.ApproximateTwoSampleKSTest
import StatsAPI.pvalue
import StatsBase.mean
import StatsBase.mad
using Plots
using DelimitedFiles
import Formatting.format

# Get relevant paths and cd to the script path.
SCRIPT_DIR = abspath(dirname(PROGRAM_FILE))
cd(SCRIPT_DIR)

function main()
    db2_path = joinpath(SCRIPT_DIR, "output", "vm2", "sweep_db_gathered.sqlite")
    if !ispath(db2_path)
        error("`output/vm2/sweep_db_gathered.sqlite` does not exist; please run gather-output.jl first")
    end
    
    db3_path = joinpath(SCRIPT_DIR, "output", "vm3", "sweep_db_gathered.sqlite")
    if !ispath(db2_path)
        error("`output/vm3/sweep_db_gathered.sqlite` does not exist; please run gather-output.jl first")
    end
    
    db2 = SQLite.DB(db2_path)
    db3 = SQLite.DB(db3_path)
    
    ts = [360 * year for year in 0:100]

    for t in ts
        println("Comparing at t = $(t)...")
        results = [
#             compare_moi([t], db2, db3, false),
#             compare_prevalence([t], db2, db3, false),
#             compare_nswitch_total([t], db2, db3, false),
#             compare_nswitch_immune([t], db2, db3, false),
#             compare_nswitch_not_immune([t], db2, db3, false),
#             compare_tswitch_sum([t], db2, db3, false),
#             compare_tswitch_overall([t], db2, db3, false),
#             compare_tswitch_not_immune([t], db2, db3, false),
#             compare_n_bites_total([t], db2, db3, false),
#             compare_n_circulating_genes([t], db2, db3, false),
#             compare_n_circulating_strains([t], db2, db3, false),
            compare_n_infected_bites([t], db2, db3, false),
#             compare_n_bites([t], db2, db3, false),
            compare_mean_age([t], db2, db3, false),
            compare_sd_age([t], db2, db3, false),
        ]
        for (title, meanmeans2, madmeans2, meanmeans3, madmeans3, pvalue) in results
            println("-----------")
            if pvalue > 0.01
                println("$(title) matches: same distribution in varmodel2 and varmodel3")
            else
                println("$(title) DOES NOT match")
            end
            println("    p-value = $(format("{:.3e}", pvalue))")
        end
        println()
    end

    ts_end = [360 * year for year in 91:100]
    results = [
        compare_moi(ts_end, db2, db3, true),
        compare_prevalence(ts_end, db2, db3, true),
        compare_nswitch_total(ts_end, db2, db3, true),
        compare_nswitch_immune(ts_end, db2, db3, true),
        compare_nswitch_not_immune(ts_end, db2, db3, true),
        compare_tswitch_sum(ts_end, db2, db3, true),
        compare_tswitch_overall(ts_end, db2, db3, true),
        compare_tswitch_not_immune(ts_end, db2, db3, true),
        compare_n_bites_total(ts_end, db2, db3, true),
        compare_n_circulating_genes(ts_end, db2, db3, true),
        compare_n_circulating_strains(ts_end, db2, db3, true),
        compare_n_infected_bites(ts_end, db2, db3, false),
#         compare_n_bites(ts_end, db2, db3, false),
        compare_mean_age(ts_end, db2, db3, false),
        compare_sd_age(ts_end, db2, db3, false),
    ]
    
    # Summary
    println("SUMMARY")
    for (title, meameans2, madmeans2, meanmeans3, madmeans3, pvalue) in results
        println("-----------")
        if pvalue > 0.01
            println("$(title) matches: same distribution in varmodel2 and varmodel3")
        else
            println("$(title) does not match")
        end
        println("    p-value = $(format("{:.3e}", pvalue))")
    end
    
    # Save p-values to CSV
    csv_data = vcat([("measurement", "meanmeans2", "madmeans2", "meanmeans3", "madmeans3", "pvalue")], results)
    writedlm("output/compare-results.csv", csv_data)
end

function compare_funcs(name, ts, db2, db3, should_plot, f2, f3)
    compare_mats(
        name, ts,
        get_matrix(ts, db2, f2),
        get_matrix(ts, db3, f3),
        should_plot
    )
end

function compare_mats(name, ts, mat2, mat3, should_plot)
    n_reps_2 = size(mat2)[2]
    n_reps_3 = size(mat3)[2]
    
    means2 = mean.(eachcol(mat2))
    meanmeans2 = mean(means2)
    madmeans2 = mad(means2; center = meanmeans2)
    means3 = mean.(eachcol(mat3))
    meanmeans3 = mean(means3)
    madmeans3= mad(means3; center = meanmeans3)
    
    if should_plot
        p = plot(ts, hcat(mat2, mat3), lcolor = reshape(vcat(
                repeat([(:black)], n_reps_2), repeat([(:green)], n_reps_3)
            ), 1, n_reps_2 + n_reps_3), lalpha = 0.3, label = false, plot_title = name
        )
        hline!(p, [meanmeans2], line = (:solid, :black), lw = 2, label = "varmodel2 ($(format("{:.2e} ({:.2e})", meanmeans2, madmeans2)))")
        hline!(p, [meanmeans3], line = (:solid, :green), lw = 2, label = "varmodel3 ($(format("{:.2e} ({:.2e})", meanmeans3, madmeans3)))")
        savefig(p, "output/$(name).pdf")
    end
    
    # Test average prevalence during years 90-100
    test = ApproximateTwoSampleKSTest(means2, means3)
#     if should_plot
#         println("Kolmogorov-Smirnov Test of $(name)...")
#         println("--------")
#         println("test results:")
#         print(test)
#         println("--------")
#         println()
#     end
    
    (name, meanmeans2, madmeans2, meanmeans3, madmeans3, pvalue(test))
end

function get_matrix(ts, db, f)
    run_ids = [run_id for (run_id,) in execute(db, "SELECT DISTINCT run_id FROM debug_stats ORDER BY run_id")]
    [f(t, db, run_id) for t in ts, run_id in run_ids]
end
    
function compare_moi(ts, db2, db3, should_plot)
    compare_funcs("moi", ts, db2, db3, should_plot,
        (
            function(t, db, run_id)
                for (n_infections, n_infected) in execute(db,
                    "SELECT n_infections, n_infected FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_infections / n_infected
                end
            end
        ),
        (
            function(t, db, run_id)
                    for (n_infections, n_infected) in execute(db3,
                        "SELECT n_infections_active, n_infected_active FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                    )
                        return n_infections / n_infected
                    end
            end
        )
    )
end

function compare_prevalence(ts, db2, db3, should_plot)
    compare_funcs("prevalence", ts, db2, db3, should_plot,
        (
            function(t, db, run_id)
                for (n_infected,) in execute(db,
                    "SELECT n_infected FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_infected / N_HOSTS
                end
            end
        ),
        (
            function(t, db, run_id)
                for (n_infected,) in execute(db3,
                    "SELECT n_infected_active FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_infected / N_HOSTS
                end
            end
        )
    )
end

function compare_nswitch_total(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (n_switch_immune, n_switch_not_immune) in execute(db, 
            "SELECT n_switch_immune, n_switch_not_immune FROM debug_stats WHERE time = ? AND run_id = ?", [t, run_id]
        )
            return Float64(n_switch_immune + n_switch_not_immune)
        end
    end
    compare_funcs("nswitch_total", ts, db2, db3, should_plot, f, f)
end

function compare_nswitch_immune(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (n_switch_immune,) in execute(db, 
            "SELECT n_switch_immune FROM debug_stats WHERE time = ? AND run_id = ?", [t, run_id]
        )
            return Float64(n_switch_immune)
        end
    end
    compare_funcs("nswitch_immune", ts, db2, db3, should_plot, f, f)
end

function compare_nswitch_not_immune(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (n_switch_not_immune,) in execute(db, 
            "SELECT n_switch_not_immune FROM debug_stats WHERE time = ? AND run_id = ?", [t, run_id]
        )
            return Float64(n_switch_not_immune)
        end
    end
    compare_funcs("n_switch_not_immune", ts, db2, db3, should_plot, f, f)
end

function compare_tswitch_sum(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (t_switch_sum,) in execute(db, 
            "SELECT t_switch_sum FROM debug_stats WHERE time = ? AND run_id = ?", [t, run_id]
        )
            return t_switch_sum
        end
    end
    compare_funcs("tswitch_sum", ts, db2, db3, should_plot, f, f)
end

function compare_tswitch_overall(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (n_switch_immune, n_switch_not_immune, t_switch_sum) in execute(db, 
            "SELECT n_switch_immune, n_switch_not_immune, t_switch_sum FROM debug_stats WHERE time = ? AND run_id = ?", [t, run_id]
        )
            return t_switch_sum / (n_switch_immune + n_switch_not_immune)
        end
    end
    compare_funcs("tswitch_overall", ts, db2, db3, should_plot, f, f)
end

function compare_tswitch_not_immune(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (n_switch_not_immune, t_switch_sum) in execute(db, 
            "SELECT n_switch_not_immune, t_switch_sum FROM debug_stats WHERE time = ? AND run_id = ?", [t, run_id]
        )
            return t_switch_sum / n_switch_not_immune
        end
    end
    compare_funcs("tswitch_not_immune", ts, db2, db3, should_plot, f, f)
end

function compare_n_bites_total(ts, db2, db3, should_plot)
    compare_funcs("n_bites_total", ts, db2, db3, should_plot,
        (
            function(t, db, run_id)
                for (n_bites,) in execute(db,
                    "SELECT n_total_bites FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_bites
                end
            end
        ),
        (
            function(t, db, run_id)
                for (n_bites,) in execute(db3,
                    "SELECT n_bites FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_bites
                end
            end
        )
    )
end

function compare_n_circulating_genes(ts, db2, db3, should_plot)
    compare_funcs("n_circulating_genes", ts, db2, db3, should_plot,
        (
            function(t, db, run_id)
                for (n_circulating_genes,) in execute(db,
                    "SELECT n_circulating_genes FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_circulating_genes
                end
            end
        ),
        (
            function(t, db, run_id)
                for (n_circulating_genes,) in execute(db3,
                    "SELECT n_circulating_genes FROM gene_strain_counts WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_circulating_genes
                end
            end
        )
    )
end

function compare_n_circulating_strains(ts, db2, db3, should_plot)
    compare_funcs("n_circulating_strains", ts, db2, db3, should_plot,
        (
            function(t, db, run_id)
                for (n_circulating_strains,) in execute(db,
                    "SELECT n_circulating_strains FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_circulating_strains
                end
            end
        ),
        (
            function(t, db, run_id)
                for (n_circulating_strains,) in execute(db3,
                    "SELECT n_circulating_strains FROM gene_strain_counts WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_circulating_strains
                end
            end
        )
    )
end

function compare_n_bites(ts, db2, db3, should_plot)
    compare_funcs("n_bites", ts, db2, db3, should_plot,
        (
            function(t, db, run_id)
                for (n_bites,) in execute(db,
                    "SELECT n_total_bites FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_bites
                end
            end
        ),
        (
            function(t, db, run_id)
                for (n_bites,) in execute(db3,
                    "SELECT n_bites FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
                )
                    return n_bites
                end
            end
        )
    )
end

function compare_n_infected_bites(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (n_infected_bites,) in execute(db,
                "SELECT n_infected_bites FROM summary WHERE time = ? AND run_id = ?", [t, run_id]
            )
            return n_infected_bites
        end
    end
    compare_funcs("n_infected_bites", ts, db2, db3, should_plot, f, f)
end

function compare_mean_age(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (mean_age,) in execute(db,
                "SELECT mean_age FROM debug_stats WHERE time = ? AND run_id = ?", [t, run_id]
            )
            return mean_age
        end
    end
    compare_funcs("mean_age", ts, db2, db3, should_plot, f, f)
end

function compare_sd_age(ts, db2, db3, should_plot)
    f = function(t, db, run_id)
        for (sd_age,) in execute(db,
                "SELECT sd_age FROM debug_stats WHERE time = ? AND run_id = ?", [t, run_id]
            )
            return sd_age
        end
    end
    compare_funcs("sd_age", ts, db2, db3, should_plot, f, f)
end

main()
