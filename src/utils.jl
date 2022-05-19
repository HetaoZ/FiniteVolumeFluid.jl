@inline function expand(v::Vector, new_length::Int, default_value) 
    return [i>length(v) ? default_value : v[i] for i in 1:new_length]
end

@inline function betweeneq(a::Vector{T}, lo::Vector{T}, hi::Vector{T}) where T <: Real
    return all(lo .≤ a .≤ hi)
end

@inline function add_cartesian(id::CartesianIndex, axis::Int, n::Int)
    return CartesianIndex(ntuple(k -> k == axis ? id[k]+n : id[k], length(id)))
end

@inline function reset_cartesian(id::CartesianIndex, axis::Int, new_index::Int)
    return CartesianIndex(ntuple(k -> k == axis ? new_index : id[k], length(id)))
end

