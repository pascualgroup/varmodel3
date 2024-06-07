using Test

import Base.push!
import Base.delete!
import Base.rand
import Random.rand!
import Random.AbstractRNG
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
function direct_sample_linear_scan(rng, weights, total_weight)
    bin_right_side = 0.0
    sampled_location = rand(rng) * total_weight
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
function shuffle_columns_to!(rng, dst, src)
    m = size(dst)[2]
    src_indices = MVector{m, Int}(1:m)
    shuffle!(rng, src_indices)
    dst[:,:] = src[:,src_indices]
end

"""
Copies a random sample of the columns of `src1` and `src2` to `dst`.
"""

function sample_columns_from_two_matrices_to_util!(rng, dst, src1, src2)
    m_dst = size(dst)[2]
    m_src_1 = size(src1)[2]
    m_src_2 = size(src2)[2]
    m_src = m_src_1 + m_src_2
    src_indices = MVector{m_src, Int}(1:m_src)
    shuffle!(rng, src_indices)
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



function sample_columns_from_two_matrices_to_util2!(rng, dst, src1, src2, P, s, infection_genes_index_var_groups)
    if !P.var_groups_fix_ratio
        sample_columns_from_two_matrices_to_util!(rng, dst, src1, src2)
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
            sample_columns_from_two_matrices_to_util!(rng, dst_genes_group_id, src1_genes_group_id, src2_genes_group_id)
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
State for a weighted discrete distribution
"""
mutable struct WeightedDiscreteDistribution
    weights::Vector{Float64}
    total_weight::Float64

    function WeightedDiscreteDistribution(bin_size, weights; rand_batch_size = 10000000)
        wdd = new(weights, sum(weights))
        wdd
    end
end

function Base.rand(rng::AbstractRNG, wdd::WeightedDiscreteDistribution)
    direct_sample_linear_scan(rng, wdd.weights, wdd.total_weight)
end

function total_weight(wdd::WeightedDiscreteDistribution)
    wdd.total_weight
end

function recompute_total_weight!(wdd::WeightedDiscreteDistribution)
    wdd.total_weight = sum(wdd.weights)
end

function update!(wdd::WeightedDiscreteDistribution, item, weight)
    if weight != wdd.weights[item]
        wdd.total_weight += weight - wdd.weights[item] # May introduce error; periodically call recompute_total_weight!()
        wdd.weights[item] = weight
    end
end

"""
"""
