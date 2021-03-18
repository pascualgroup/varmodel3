#!/usr/bin/env julia

"""
This script performs a single run of model from a parameters file in JSON format.

To run it, provide the location of the parameters file as the first command-line
argument. (If not provided, "parameters.json" is assumed.)

```
<path1>/run.jl <path2>/parameters.json
```
"""

using JSON3

include("preamble.jl")

# Load parameters into a constant global, so that they can be used when compiling
# the model code.
const P = let
    params_filename = if length(ARGS) == 0
        "parameters.json"
    elseif length(ARGS) == 1
        ARGS[1]
    else
        error("Usage: <path-to>/run.jl [parameters.json]")
    end
    
    json_str = read(params_filename, String)
    JSON3.read(json_str, Params)
end

include("src/model.jl")
run()
