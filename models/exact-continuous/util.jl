function direct_sample_linear_scan(weights, total_weight)
    bin_right_side = 0.0
    sampled_location = rand() * total_weight
    for i in 1:(length(weights) - 1)
        bin_right_side += weights[i]
        if sampled_location < bin_right_side
            return i
        end
    end
    length(weights)
end
