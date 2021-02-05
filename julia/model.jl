include("util.jl")
include("exact.jl")

function run(p::Params)
    validate(p)
    
    if p.implementation == EXACT_GILLESPIE
        run_exact(p)
    else
        @assert false
    end
end

