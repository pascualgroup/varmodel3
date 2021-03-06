using Test

import Base.push!
import Base.delete!
import Base.rand
import Base.length
import Base.iterate
import Base.in

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
