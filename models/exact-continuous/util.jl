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

function shuffle_to!(dst, src)
    copyto!(dst, src)
    shuffle!(dst)
    nothing
end

function sample_from_two_vectors_to!(dst, src1, src2)
    m_dst = length(dst)
    m_src_1 = length(src1)
    m_src_2 = length(src2)
    m_src = m_src_1 + m_src_2
    src_indices = MVector{m_src, Int}(1:m_src)
    shuffle!(src_indices)
    for i_dst in 1:m_dst
        i_src = src_indices[i_dst]
        if i_src <= m_src_1
            dst[i_dst] = src1[i_src]
        else
            dst[i_dst] = src2[i_src - m_src_1]
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
