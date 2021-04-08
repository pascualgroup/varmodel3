#!/usr/bin/env julia

"""
This script runs multiple external commands in parallel.

Usage:

```
runmany.jl <n-cores> <command-list-filename>
```

The program reads a list of commands from the lines in
`<command-list-filename>` and executes them using Julia's built-in `asyncmap()`
function using `<n-cores>` tasks. Because each command is spawned in a separate
process, `<n-cores>` of them will execute in parallel at any given time.
As commands finish, remaining commands will be dispatched from the list.

All commands are executed in the working directory this script was executed in.
If they need to operate in a different directory, they need to change to that
directory themselves.
"""

function main()
    n_cores = parse(Int, ARGS[1])
    run_cmds_filename = ARGS[2]
    
    run_cmds = readlines(run_cmds_filename)
    results = asyncmap(do_run, run_cmds; ntasks = n_cores)
    nothing
end
    
function do_run(run_cmd)
    println("$(run_cmd)\n  starting")
    run(Cmd(`$(run_cmd)`; ignorestatus = true))
    println("$(run_cmd)\n  done")
end

main()
