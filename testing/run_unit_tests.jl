#!/usr/bin/env julia

using Test

@testset "Unit tests for varmodel3" begin
    @testset "Utility tests" begin
        include("util_tests.jl")
    end
end