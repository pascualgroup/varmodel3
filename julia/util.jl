import Base.push!
import Base.delete!
import Base.rand
import Base.length
import Base.iterate

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

function push!(c::IndexedSet{T}, x::T) where {T}
    push!(c.array, x)
    c.indices[x] = lastindex(c.array)
end

function delete!(c::IndexedSet{T}, x::T) where {T}
    if !haskey(c.indices, x)
        return c
    end
    
    a = c.array
    index = c.indices[x]
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
