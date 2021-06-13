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
