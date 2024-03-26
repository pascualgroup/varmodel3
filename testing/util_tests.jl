#!/usr/bin/env julia

using Test
import HypothesisTests.ChisqTest
import HypothesisTests.pvalue
import Distributions.LogNormal

include("../src/util.jl")

@testset "WeightedDiscreteDistribution tests" begin
    function pvalue_wdd(wdd, n_draws)
        counts = zeros(Int, length(wdd.weights))
        for i in 1:n_draws
            counts[rand(wdd)] += 1
        end
        println(wdd.weights)
        println(counts)
        probs = wdd.weights / sum(wdd.weights)
        println(probs)
        test = ChisqTest(counts, probs)
        println(test)
        pvalue(test)
    end

    function pvalue_new_wdd(bin_size, weights, n_draws)
        wdd = WeightedDiscreteDistribution(bin_size, weights)
        pvalue_wdd(wdd, n_draws)
    end

    function fuzz_modify_weights(wdd, weight_dist, n_iterations)
        for i in 1:n_iterations
            item = rand(1:length(wdd.weights))
            update!(wdd, item, rand(weight_dist))
        end
    end

    @testset "Static tests" begin
        @test pvalue_new_wdd(0.1, [0.29, 1.23, 127.0], 100000) > 0.01
        @test pvalue_new_wdd(0.1, [0.1, 0.2, 0.4, 1.0], 100000) > 0.01
        @test pvalue_new_wdd(100.0, [0.1, 0.2, 99.0], 100000) > 0.01
    end

    @testset "Dynamic tests" begin
        @testset "Dynamic Test 1" begin
            weight_dist = LogNormal(1.0, 1.0)
            weights_initial = rand(weight_dist, 4)
            wdd = WeightedDiscreteDistribution(0.1, weights_initial)
            @test pvalue_wdd(wdd, 100000) > 0.01
            fuzz_modify_weights(wdd, weight_dist, 10000)
            @test pvalue_wdd(wdd, 100000) > 0.01
        end
    end
end
