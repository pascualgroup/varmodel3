using Test

import Base.push!
import Base.delete!
import Base.rand
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
function sample_columns_from_two_matrices_to!(dst, src1, src2)
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

"""
Deletes index `i` in the array and replaces it with the last entry in the array.

This enables maintaining a set of objects in an array (in arbitrary order)
so they can be uniformly randomly sampled, with constant-time cost for removal.
"""
function delete_and_swap_with_end!(a, i)
    @assert i <= length(a)

    x = pop!(a)
    if i <= length(a)
        a[i] = x
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
