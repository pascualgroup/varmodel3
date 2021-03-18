# Check that "parameters.jl" got loaded already, and that the global parameters
# constant P is defined
@assert @isdefined Params
@assert @isdefined P
@assert typeof(P) === Params
validate(P)

# Load shared files
include("output.jl")
include("util.jl")

# Load code for the model variant/implementation specified in the parameters
if P.use_discrete_time_approximation
    # Discrete-time code needs some updating to be consistent with the reorganization
    @assert false
else
    include("continuous/model.jl")
end
