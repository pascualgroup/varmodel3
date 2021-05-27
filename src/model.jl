"""
This file loads code shared by all model implementations, and then loads the
variant of the model specified by the global parameters constant P.

The `run()` function, called in run scripts, is defined by the specific variant
of the model that is loaded, e.g., inside `continuous/model.jl`.

This is an unsophisticated but straightforward way to have multiple variants of
the model share a parameters format, run script, and other bits of code.
"""

# Verify that "parameters.jl" got loaded already, and that the global parameters
# constant P is defined
@assert @isdefined Params
@assert @isdefined P
@assert typeof(P) === Params
validate(P)

# Load shared files
include("output.jl")
include("util.jl")

# Load code for the model specified in the parameters
if P.model == VAR_ONLY
    include("var_only/model.jl")
else
    include("var_with_parasitemia/model.jl")
end
