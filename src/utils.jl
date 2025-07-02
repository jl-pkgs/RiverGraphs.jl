export unlist, reverse_index
# import Ipaper: unlist

unlist(x::Real) = x
function unlist(p::Pair)
  [map(unlist, p.first) |> x -> vcat(x...); p.second] |> sort
end

# min_sto = get_MinSto(streamorder; level, min_sto)
function get_MinSto(streamorder; level::Int=2, min_sto=nothing)
  isnothing(min_sto) && (min_sto = max(maximum(streamorder) - level, 1))
  return min_sto
end


"""
    active_indices(subcatch_2d, nodata)

Takes a 2D array of the subcatchments. And derive forward and reverse indices.

1: Get a list of `CartesianIndex{2}`` that are active, based on a nodata value.
These map from the 1D internal domain to the 2D external domain.

2: Make a reverse index, a `Matrix{Int}``, which maps from the 2D external domain to
the 1D internal domain, providing an Int which can be used as a linear index. Values of 0
represent inactive cells.
"""
function active_indices(A::AbstractMatrix, nodata)
  all_inds = CartesianIndices(size(A))
  indices = filter(i -> !isequal(A[i], nodata), all_inds) # 不为NA的全部收纳

  reverse_indices = zeros(Int, size(A))
  for (i, I) in enumerate(indices)
    reverse_indices[I] = i
  end
  return indices, reverse_indices
end


function reverse_index(v::AbstractVector{T}, index, index_rev; nodata::T=T(0)) where {T}
  R = zeros(T, size(index_rev)) .+ nodata
  @inbounds for (i, I) in enumerate(index)
    R[I] = v[i]
  end
  R
end

function Base.Matrix(g::RiverGraph, v::AbstractVector{T}, nodata::T=T(0)) where {T}
  reverse_index(v, g.index, g.index_rev; nodata)
end

function Base.Matrix(g::RiverGraph, xs::Tuple, nodata=0)
  map(v -> Matrix(g, v, nodata), xs)
end

Base.Matrix(g::RiverGraph) = Matrix(g, g.data, g.nodata)
