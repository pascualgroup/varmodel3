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

function get_key_by_iteration_order(d::Dict{K, V}, index::Int) where {K, V}
    for (i, k) in enumerate(keys(d))
        if i == index
            return k
        end
    end
    @assert false
end
