include("parameters.jl")
include("output.jl")
include("util.jl")
include("shared.jl")
include("discrete_time.jl")

function run(p::Params)
    validate(p)
    
    if p.use_discrete_time_approximation
        run_discrete_time(p)
    else
        @assert false
    end
end
