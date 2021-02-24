include("util.jl")
include("shared.jl")
include("discrete_time.jl")

function run(p::Params)
    validate(p)
    
    if p.implementation == CONTINUOUS_TIME
        # run_continuous_time(p)
    elseif p.implementation == DISCRETE_TIME
        run_discrete_time(p)
    else
        @assert false
    end
end
