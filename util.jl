using Test

import Base.push!
import Base.delete!
import Base.rand
import Base.length
import Base.iterate
import Base.in

using Distributions

function test_utils()
    @testset begin
        test_mask_with_row_limits()
    end
end

"""
    A random-samplable mutable set.
    
    Stored as an array, which is in arbitrary order, and a dictionary of indices,
    which are used to efficiently remove items from the set.
"""
struct IndexedSet{T}
    array::Array{T}
    indices::Dict{T, Int}
    
    function IndexedSet{T}() where {T}
        new(T[], Dict{T, Int}())
    end
end

Base.eltype(::Type{IndexedSet{T}}) where {T} = T

function in(x::T, c::IndexedSet{T}) where {T}
    haskey(c.indices, x)
end

function push!(c::IndexedSet{T}, x::T) where {T}
    push!(c.array, x)
    c.indices[x] = lastindex(c.array)
end

function delete!(c::IndexedSet{T}, x::T) where {T}
#     if !haskey(c.indices, x)
#         return c
#     end
    
    a = c.array
    index = c.indices[x]
    delete!(c.indices, x)
    if index != lastindex(a)
        y = last(a)
        a[index] = y
        c.indices[y] = index
    end
    pop!(a)
    c
end

function length(c::IndexedSet{T}) where {T}
    length(c.array)
end

function rand(c::IndexedSet{T}) where {T}
    rand(c.array)
end

function delete!(a::Vector, x)
    for i in 1:lastindex(a)
        if a[i] == x
            return deleteat!(a, i)
        end
    end
end

function iterate(c::IndexedSet{T}, state = 1) where {T}
    state > length(c) ? nothing : ( c.array[state], state + 1 )
end

function findall!(x::Array{Int}, y)
    empty!(x)
    for (i, yi) in enumerate(y)
        if yi
            push!(x, i)
        end
    end
end

function collect!(x, y)
    empty!(x)
    append!(x, y)
end

function shuffleto!(x, y)
    collect!(x, y)
    shuffle!(x)
end

function vcat!(x, y, z)
    empty!(x)
    append!(x, y)
    append!(x, z)
end

function mask_with_row_limits(mask, limits)
#     println("mask = $(mask), limits = $(limits)")
    @assert length(limits) == size(mask)[1]
    count = fill(0, length(limits))
    new_mask = falses(size(mask))
    for i in 1:size(mask)[2]
        new_mask[:,i] = mask[:,i] .& (count .< limits)
        count[:] .+= new_mask[:,i]
    end
#     println("new_mask = $(new_mask)")
    new_mask
end

function test_mask_with_row_limits()
    testval = mask_with_row_limits(trues(2, 1), [1, 0]) == Bool[1; 0]
    @testset begin
        @test begin
            mask_with_row_limits(falses(1, 1), [0]) == falses(1, 1)
        end
        @test begin
            mask_with_row_limits(trues(1, 1), [0]) == falses(1, 1)
        end
        @test begin
            mask_with_row_limits(trues(1, 1), [1]) == trues(1, 1)
        end
        @test begin
            mask_with_row_limits(falses(2, 1), [0, 0]) == falses(2, 1)
        end
        @test begin
            mask_with_row_limits(trues(2, 1), [0, 0]) == falses(2, 1)
        end
        @test begin
            mask_with_row_limits(trues(2, 1), [1, 0]) == reshape([true false], (2, 1))
        end
        @test begin
            mask_with_row_limits(trues(2, 2), [1, 0]) == [true false; false false]
        end
    end
end

function fill_mask_with_entries(mask, entries)
    x = fill(entries[1], size(mask))
    x[mask] .= entries
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

function sample_each_column_without_replacement(x, n)
    samples = fill(x[1], (n, size(x)[2]))
    for j in 1:size(x)[2]
        samples[:,j] = sample(x[:,j], n; replace = false)
    end
    samples
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

function zip_index(x, i...)
    x[zip_cartesian(i...)]
end

function zip_view(x, i...)
    @view x[zip_cartesian(i...)]
end

"""
    Identify columns already present in a matrix, and return column indices
    for each.
    
    Zeros in the returned indices indicate missing columns.
"""
function match_columns(mat, cols)
#     println("size(mat) = $(size(mat))")
#     println("size(cols) = $(size(cols))")
    indices = fill(0, size(cols)[2])
    for j in 1:size(cols)[2]
        index = findfirst(reshape(
            all(mat .== cols[:, j:j]; dims = 1),
            size(mat)[2]
        ))
        if index != nothing
            indices[j] = index
        end
    end
    indices
end

"""
    Remove non-monotonic increases
"""
function remove_index_gaps(x)
    xp = fill(0, size(x))
    index_map = fill(0, size(x))
    last_max_x = 0
    last_max_xp = 0
    for (i, j) in enumerate(x)
        if j > last_max_x
            last_max_x = j
            last_max_xp += 1
            index_map[last_max_x] = last_max_xp
        end
        xp[i] = index_map[j]
    end
    xp
end

"""
    An array of dictionaries stored using large arrays with stable memory
    layout. Keys are length-D tuples of positive integers of type K,
    and values are of type V.
    
    The dimensions represent:
    (1) Entries in the hash table corresponding to the same bucket (index).
    (2) Indices in the hash table corresponding to (hash value) % (size)
    (3) Index in the array of dictionaries.
    
    The structure is dynamically resized by 25% in dimension 1 if the number of
    entries exceeds the current capacity for any single bucket.
    
    The structure is dynamically resized by 25% in dimension 2 if the load
    factor exceeds 0.75. This triggers reallocation and rehashing of all entries.
    The load factor is computed as an average across all dicts in the array.
"""
mutable struct ArrayOfTupleDicts{K, D, V}
    keys::Array{K, 4}
    values::Array{V, 3}
    
    function ArrayOfTupleDicts{K, D, V}(n_dicts, n_entries_per_dict_guess) where {K, D, V}
        n = n_entries_per_dict_guess
        k = (n + 1) * 4 ÷ 3
        max_n_per_k = guess_max_entries_per_bucket(n, k)
        
        new{K, D, V}(
            fill(K(0), (D, max_n_per_k, k, n_dicts)),
            fill(V(0), (max_n_per_k, k, n_dicts))
        )
    end
end

function dict_count(x::ArrayOfTupleDicts{K, D, V}) where {K, D, V}
    size(x.values)[3]
end

function bucket_count(x::ArrayOfTupleDicts{K, D, V}) where {K, D, V}
    size(x.values)[2]
end

function entries_per_bucket(x::ArrayOfTupleDicts{K, D, V}) where {K, D, V}
    size(x.values)[1]
end

function tuple_size(x::ArrayOfTupleDicts{K, D, V}) where {K, D, V}
    size(x.keys)[1]
end

function guess_max_entries_per_bucket(n, k)
    # The distribution of entries in a single bucket is Binomial(n, 1/k).
    # We'll guess max entries per bucket as 25% more than the number of entries that
    # reaches the 1 / n quantile of the CDF of this distribution.
    # 
    # The correct guess would be, e.g., the median of the distribution of the
    # *maximum* number of entries across all n entries, but that seems
    # hard to figure out.
    
    distribution = Binomial(n, 1.0 / k)
    target_quantile = (n - 1) / n
    entries_per_bucket = 1
    while true
        if cdf(distribution, entries_per_bucket) > target_quantile
            break
        end
        entries_per_bucket *= 2
    end
    (entries_per_bucket + 1) * 5 ÷ 4
end


