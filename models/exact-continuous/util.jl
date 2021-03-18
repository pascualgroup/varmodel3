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
    copyto!(dst, src)
    shuffle_columns!(dst)
    nothing
end

function shuffle_columns!(mat)
    m = size(mat)[2]
    for i in 1:(m - 1)
        tmp = mat[:, i]
        j = rand((i+1):m)
        mat[:, i] = mat[:, j]
        mat[:, j] = tmp
    end
    nothing
end

# function xor_swap!(x, y)
#     x .⊻= y
#     y .⊻= x
#     x .⊻= y
#     nothing
# end

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