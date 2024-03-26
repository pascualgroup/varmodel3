using Test

import Base.push!
import Base.delete!
import Base.rand
import Random.rand!
import Base.length
import Base.iterate
import Base.in

using Distributions

"""
Samples an index using weights.

Assumes that `total_weight == sum(weights)`.

The sampling algorithm can be visualized by arranging boxes left to right,
with widths equal to the weights, in sequence, numbered `1:length(weights)`;
the leftmost box has its left edge at 0.

```
|----------------------------------------------------|
| weights[1] |     weights[2]     |    weights[3]    |
|-----------------|----------------------------------|
0           sampled_location                    total_weight
```

To sample, draw a real number is uniformly randomly in [0, total_weight];
this sampled location corresponds to a point within one of the boxes.
To identify which box, scan the boxes left to right while computing the
cumulative sum.
If the sampled location is less than the cumulative sum, return the index of
the most recently scanned box.
"""
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

"""
Copies the columns of `src` to `dst` in random order.
"""
function shuffle_columns_to!(dst, src)
    m = size(dst)[2]
    src_indices = MVector{m, Int}(1:m)
    shuffle!(src_indices)
    dst[:,:] = src[:,src_indices]
end

"""
Copies a random sample of the columns of `src1` and `src2` to `dst`.
"""

function sample_columns_from_two_matrices_to_util!(dst, src1, src2)
    m_dst = size(dst)[2]
    m_src_1 = size(src1)[2]
    m_src_2 = size(src2)[2]
    m_src = m_src_1 + m_src_2
    src_indices = MVector{m_src, Int}(1:m_src)
    shuffle!(src_indices)
    for i_dst in 1:m_dst
        i_src = src_indices[i_dst]
        if i_src <= m_src_1
            dst[:, i_dst] = src1[:, i_src]
        else
            dst[:, i_dst] = src2[:, i_src - m_src_1]
        end
    end
    nothing
end



function sample_columns_from_two_matrices_to_util2!(dst, src1, src2, P, s, infection_genes_index_var_groups)
    if !P.var_groups_fix_ratio
        sample_columns_from_two_matrices_to_util!(dst, src1, src2)
    else
        # extract genes of the specific group from two source strains.
        src1_genes_all_groups = []
        src2_genes_all_groups = []
        for group_id in 1:length(P.var_groups_ratio)
            push!(src1_genes_all_groups, [])
            push!(src2_genes_all_groups, [])
        end
        
        for i in 1:size(src1)[2]
            src1_gene = Gene(src1[:, i])
            # @assert haskey(s.association_genes_to_var_groups, src1_gene)
            src1_gene_group_id = s.association_genes_to_var_groups[src1_gene]
            push!(src1_genes_all_groups[src1_gene_group_id], src1[:, i])
        end
        for i in 1:size(src2)[2]    
            src2_gene = Gene(src2[:, i])
            # @assert haskey(s.association_genes_to_var_groups, src2_gene)
            src2_gene_group_id = s.association_genes_to_var_groups[src2_gene]
            push!(src2_genes_all_groups[src2_gene_group_id], src2[:, i])
        end
        
        for group_id in 1:length(P.var_groups_ratio)
            dst_index = infection_genes_index_var_groups[group_id]
            src1_genes_group_id_vec = src1_genes_all_groups[group_id]
            # @assert length(src1_genes_group_id_vec) == round(Int, P.var_groups_ratio[group_id] * P.n_genes_per_strain) 
            src1_genes_group_id = zeros(AlleleId, P.n_loci, length(src1_genes_group_id_vec))
            for j in 1:length(src1_genes_group_id_vec)
                src1_genes_group_id[:,j] = src1_genes_group_id_vec[j]
            end
            src2_genes_group_id_vec = src2_genes_all_groups[group_id]
            # @assert length(src2_genes_group_id_vec) == round(Int, P.var_groups_ratio[group_id] * P.n_genes_per_strain) 
            src2_genes_group_id = zeros(AlleleId, P.n_loci, length(src2_genes_group_id_vec))
            for k in 1:length(src2_genes_group_id_vec)
                src2_genes_group_id[:,k] = src2_genes_group_id_vec[k]
            end
            dst_genes_group_id = zeros(AlleleId, P.n_loci, length(dst_index))
            sample_columns_from_two_matrices_to_util!(dst_genes_group_id, src1_genes_group_id, src2_genes_group_id)
            dst[:, dst_index] = dst_genes_group_id
        end
    end
    nothing
end

"""
Calculates the duration of an infection.
"""
function get_duration!(a, i, t)
    @assert i <= length(a)
    if i <= length(a)
        a[i].duration = t - a[i].t_infection
    end
end

"""
Deletes index `i` in the array and replaces it with the last entry in the array.

This enables maintaining a set of objects in an array (in arbitrary order)
so they can be uniformly randomly sampled, with constant-time cost for removal.
The moved item is returned so an auxiliary data structure can also track locations if necessary.
"""
function delete_and_swap_with_end!(a, i)
    @assert i <= length(a)
    
    x = pop!(a)
    if i <= length(a)
        a[i] = x
        return x
    end
    nothing
end

"""
Returns a key in a dictionary, by iteration order.

Dictionary iteration order is undefined, but if `index` has been uniformly
randomly sampled, this will have the effect of uniformly randomly sampling
a key in the dictionary.

This is used to avoid the extra memory expense of maintaining a randomly
samplable array.
"""
function get_key_by_iteration_order(d::Dict{K, V}, index::Int) where {K, V}
    for (i, k) in enumerate(keys(d))
        if i == index
            return k
        end
    end
    @assert false
end

"""
State for a batched distribution.
"""
mutable struct BatchedDistribution
    d::Sampleable{Univariate, Continuous}
    draws::Vector{Float64}
    next_draw_index::Int

    function BatchedDistribution(d, batch_size)
        draws = Vector(undef, batch_size)
        next_draw_index = 1
        new(d, draws, 1)
    end
end

function Base.rand(bd::BatchedDistribution)
    if bd.next_draw_index == 1
        rand!(bd.d, bd.draws)
    end
    draw = bd.draws[bd.next_draw_index]
    bd.next_draw_index = 1 + mod(bd.next_draw_index, length(bd.draws))
    draw
end

"""
State for a weighted discrete distribution
"""
mutable struct WeightedDiscreteDistribution
    bin_size::Float64
    weights::Vector{Float64}
    bins::Vector{Int}
    bin_indices::Dict{Int, Set{Int}}
    p_accept::Vector{Float64}
    batch_dist::BatchedDistribution

    function WeightedDiscreteDistribution(bin_size, weights; rand_batch_size = 10000000)
        # Initialize bins
        bins = Vector{Int}()
        bin_indices = Dict{Int, Set{Int}}()
        p_accept = Vector{Float64}()
        for i in 1:length(weights)
            (n_bins, p_accept_i) = compute_n_bins_and_p_accept(weights[i], bin_size)
            push!(p_accept, p_accept_i)

            bin_indices[i] = Set((length(bins) + 1):(length(bins) + n_bins))
            @assert length(bin_indices[i]) == n_bins
            append!(bins, fill(i, n_bins))
        end
        #println(p_accept)

        wdd = new(bin_size, weights, bins, bin_indices, p_accept, BatchedDistribution(Uniform(), rand_batch_size))
        verify(wdd)
        wdd
    end
end

function Base.rand(wdd::WeightedDiscreteDistribution)
    # Repeat until a sample is accepted
    while true
        item = rand(wdd.bins) # TODO: batch this too if not fast enough
        if rand(wdd.batch_dist) < wdd.p_accept[item]
            return item
        end
    end
end

function compute_n_bins_and_p_accept(weight, bin_size)
    n_bins_fractional = weight / bin_size
    n_bins_ceil = ceil(n_bins_fractional)

    if n_bins_ceil == 0.0
        (0, 0.0)
    else
        (Int(n_bins_ceil), n_bins_fractional / n_bins_ceil)
    end
end

function verify(wdd::WeightedDiscreteDistribution)
    for (item, weight) in enumerate(wdd.weights)
        (n_bins, p_accept) = compute_n_bins_and_p_accept(weight, wdd.bin_size)
        @assert wdd.p_accept[item] == p_accept
        @assert length(wdd.bin_indices[item]) == n_bins
        for bin_index in wdd.bin_indices[item]
            @assert wdd.bins[bin_index] == item
        end
    end

    for (bin_index, item) in enumerate(wdd.bins)
        @assert bin_index in wdd.bin_indices[item]
    end
end

function update!(wdd::WeightedDiscreteDistribution, item, weight)
    if weight != wdd.weights[item]
        n_bins_old = length(wdd.bin_indices[item])
        (n_bins_new, p_accept_new) = compute_n_bins_and_p_accept(weight, wdd.bin_size)

        # Adjust the number of bins, if necessary
        if n_bins_new > n_bins_old
            for i in 1:(n_bins_new - n_bins_old)
                push!(wdd.bins, item)
                @assert !(length(wdd.bins) in wdd.bin_indices[item])
                push!(wdd.bin_indices[item], length(wdd.bins))
            end
        elseif n_bins_old > n_bins_new
            for i in 1:(n_bins_old - n_bins_new)
                bin_index_to_remove = pop!(wdd.bin_indices[item])
                last_bin_index = length(wdd.bins)
                moved_item = delete_and_swap_with_end!(wdd.bins, bin_index_to_remove)
                if moved_item !== nothing
                    # Update the index set for the item that got moved from last_index
                    # to bin_index_to_remove
                    @assert last_bin_index in wdd.bin_indices[moved_item]
                    delete!(wdd.bin_indices[moved_item], last_bin_index)
                    @assert !(bin_index_to_remove in wdd.bin_indices[moved_item])
                    push!(wdd.bin_indices[moved_item], bin_index_to_remove)

                end
            end
        end

        # Update the weight and acceptance probability
        wdd.weights[item] = weight
        wdd.p_accept[item] = p_accept_new

        # verify(wdd)
    end
end
