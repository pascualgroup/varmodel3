include("util.jl")
include("exact.jl")
include("discrete.jl")

function run(p::Params)
    validate(p)
    
    if p.implementation == EXACT_GILLESPIE
        run_exact(p)
    elseif p.implementation == DISCRETE_APPROXIMATION
        run_discrete(p)
    else
        @assert false
    end
end

