using Test

import Base.push!
import Base.delete!
import Base.rand
import Base.length
import Base.iterate
import Base.in

using Distributions

function fill_mask_with_entries(mask, entries::Union{AbstractArray{T}, AbstractRange{T}}) where {T}
    x = fill(T(0), size(mask))
    x[mask] = entries
    x
end

function mask_with_column_limits(mask, limits)
#     println("mask = $(mask), limits = $(limits)")
    @assert length(limits) == size(mask)[2]
    count = fill(0, length(limits))
    new_mask = falses(size(mask))
    for i in 1:size(mask)[1]
        new_mask[i,:] = mask[i,:] .& (count .< limits)
        count[:] .+= new_mask[i,:]
    end
#     println("new_mask = $(new_mask)")
    new_mask
end

function sample_true_indices_by_column(mask)
    @assert all(sum(mask; dims = 1) .> 0)
    
    n, m = size(mask)
    row_indices = fill(0, m)
    col_indices = collect(1:m)
    while length(col_indices) > 0
        new_row_indices = rand(1:n, length(col_indices))
        valid_row_indices = mask[zip_cartesian(new_row_indices, col_indices)]
        row_indices[col_indices] .= new_row_indices .* valid_row_indices
        col_indices = col_indices[.!valid_row_indices]
    end
    
    row_indices
end

function findfirst_each_column(x)
    inds = fill(0, size(x)[2])
    for j in 1:size(x)[2]
        inds[j] = findfirst(x[:,j])
    end
    inds
end

function zip_cartesian(x...)
    collect(CartesianIndex(x) for x in zip(x...))
end

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

function shuffle_columns_to!(dst, src)
    m = size(dst)[2]
    src_indices = MVector{m, Int}(1:m)
    shuffle!(src_indices)
    dst[:,:] = src[:,src_indices]
end

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

function delete_and_swap_with_end!(a, i)
    @assert i <= length(a)
    
    x = pop!(a)
    if i <= length(a)
        a[i] = x
    end
    nothing
end

function get_key_by_iteration_order(d::Dict{K, V}, index::Int) where {K, V}
    for (i, k) in enumerate(keys(d))
        if i == index
            return k
        end
    end
    @assert false
end

"""
    Interpolates y(x) from vector of y's.
    
    Assumes x-values corresponding to vector of y's are at
    (0.0 + offset, 1.0 + offset, ...)
    
    Does not extrapolate beyond endpoints.
"""
function interpolate_linear(x, ys; offset = 0.5)
    x_offset = x - offset
    x_offset_floor = floor(x_offset)
    x_offset_diff = x_offset - x_offset_floor
    index = Int(x_offset_floor) + 1
    
    if x_offset_diff == 0.0
        @assert 1 <= index <= length(ys)
        return ys[index]
    end
    
    @assert 1 <= index < length(ys)
    
    ys[index] * (1.0 - x_offset_diff) + ys[index + 1] * x_offset_diff
end

"""
    Interpolates y(x) from parallel vectors of x's and y's.
    
    Does not extrapolate beyond endpoints.
"""
function interpolate_linear(x, xs, ys; offset = 0.5)
    @assert length(xs) == length(ys)
    
    # Linear search. Switch to binary if this becomes a bottleneck.
    index = 0
    for i in 1:(length(xs) - 1)
        if xs[i] <= x < xs[i + 1]
            index = i
            break
        end
    end
    
    if index == 0 && x == xs[length(xs)]
        index = length(xs)
    end
    
    x_diff = x - xs[index]
    
    if x_diff == 0.0
        @assert 1 <= index <= length(xs)
        return ys[index]
    end
    
    @assert 1 <= index < length(ys)
    
    w = x_diff / (xs[index + 1] - xs[index])
    ys[index] * (1.0 - w) + ys[index + 1] * w
end

"""
    Append to a key vector and a value vector where keys are sorted.
    If the number of entries is `large_size` or bigger and a power of two,
    remove all the entries where `key < threshold`.
    
    This function will maintain the vectors in a state where any
    `key >= threshold` will be kept, but some `key < threshold` may be thrown
    away.
    
    This is used to ensure that we always have parasitemia history values up to
    n days ago.
"""
function push_kv_keep_ge!(keys, values, k, v, threshold; large_size = 8)
    @assert length(keys) == length(values)
    n = length(keys)
    
    # for n >= 1, (n & (n - 1) == 0) checks if n is a power of two
    if n >= large_size && (n & (n - 1) == 0)
        start = findfirst_ge(keys, threshold)
        if start != 0
            stop = length(keys)
            n_new = stop - start + 1
            
            if n_new > 0
                keys[1:n_new] = keys[start:stop]
                values[1:n_new] = values[start:stop]
            end
            
            resize!(keys, n_new)
            resize!(values, n_new)
        end
    end
    
    push!(keys, k)
    push!(values, v)
end

"""
    Find the first index in v greater than or equal to x.
"""
function findfirst_ge(v, x)
    # Linear search. Switch to binary if this becomes a bottleneck
    for (i, vi) in enumerate(v)
        if vi >= x
            return i
        end
    end
    return 0
end
