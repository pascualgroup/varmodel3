#!/usr/bin/env julia

"""
This script runs commands in parallel using worker threads.

Usage:

```
runmany.jl <n-cores> <list-of-runs>
```
"""

function main()
    n_cores = parse(Int, ARGS[1])
    run_cmds_filename = ARGS[2]
    
    run_cmds = readlines(run_cmds_filename)
    results = asyncmap(do_run, run_cmds; ntasks = n_cores)
end
    
function do_run(run_cmd)
    println("$(run_cmd)\n  starting")
    run(`$(run_cmd)`)
    println("$(run_cmd)\n  done")
end

main()
