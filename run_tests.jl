#!/usr/bin/env julia

include("parameters.jl")
include("output.jl")
include("util.jl")
include("shared.jl")
include("discrete_time.jl")

function main()
    @testset begin
        test_utils()
    end
end

main()
