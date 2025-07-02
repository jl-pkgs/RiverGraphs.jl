"""
    stream_link(g, streamorder, toposort, min_sto)

Return stream_link with a unique id starting at 1, based on a minimum streamorder `min_sto`,
directed acyclic graph `g` and topological order `toposort`.
"""
function stream_link(g, toposort, streamorder; level::Int=2, min_sto=nothing)
  min_sto = get_MinSto(streamorder; level, min_sto)
  n = length(toposort)
  links = fill(0, n)

  i = 1
  for v in toposort
    streamorder[v] < min_sto && continue

    ds_nodes = outneighbors(g, v)
    if !isempty(ds_nodes)
      if streamorder[v] != streamorder[only(ds_nodes)]
        links[v] = i
        i += 1
      end
    else
      # also set pits (without a downstream node)
      links[v] = i
      i += 1
    end
  end
  return links
end

stream_link(g::RiverGraph, strord; level::Int=2, min_sto=nothing) =
  stream_link(g.graph, g.toposort, strord; level, min_sto)



function getInfo_links(g::RiverGraph, links_2d)
  con = links_2d .!= 0
  xs, ys = get_coord(con)
  vals = filter(x -> x != 0, links_2d)
  inds = g.index_rev[con][:]
  DataFrame(; x=xs, y=ys,
    lon=g.lon[xs], lat=g.lat[ys],
    link=vals, index=inds)
end



link2index(links::AbstractVector) = findall(links .!== 0)

function link2point(rg::RiverGraph, links::AbstractVector)
  inds = link2index(links)
  index2point(rg, inds)
end

function index2link(rg::RiverGraph, index::AbstractVector)
  n = length(rg.toposort)
  links = fill(0, n)
  for (i, index) in enumerate(index)
    links[index] = i
  end
  links
end

function index2point(rg::RiverGraph, inds::AbstractVector)
  (; lon, lat) = rg
  map(I -> begin
      index = rg.index[I]
      _i = index[1]
      _j = index[2]
      (lon[_i], lat[_j])
    end, inds)
end


# point2index(rg::RiverGraph, points) = find_pits(rg, points)
# 根据经纬度，找到pits
function point2index(rg::RiverGraph, points)
  (; lon, lat) = rg
  locs = st_location_exact2(lon, lat, points) # 查找位置
  map(p -> rg.index_rev[p[1], p[2]], locs) # index_pit, 流域出水口的位置
end

function move2next(rg::RiverGraph, vs::AbstractVector{Int})
  map(v -> outneighbors(rg.graph, v) |> only, vs)
end

function move2next(rg::RiverGraph, points::Vector{Tuple{Float64,Float64}})
  index_pour = point2index(rg, points)  # 下标
  index_pit = move2next(rg, index_pour) # 由于pour flowdir not na, 往下移动一个网格
  index2point(rg, index_pit)
end


function find_pits(rg::RiverGraph, points)
  index_pour = point2index(rg, points)  # 下标
  move2next(rg, index_pour) # 由于pour flowdir not na, 往下移动一个网格
end



# 有些水文站不在河流的交汇处，需要手动加上
function add_links!(links, index_pit)
  nodes_hit = findall(!isequal(0), links)
  nodes_miss = setdiff(index_pit, nodes_hit)
  n = maximum(links)
  if !isempty(nodes_miss)
    for (i, v) = enumerate(nodes_miss)
      links[v] = n + i
    end
  end
end


export find_pits, move2next, add_links!
export stream_link, getInfo_links
export point2index,
  index2link, index2point,
  link2index, link2point



function findnear(x::Real, vals::AbstractVector; cell::Real=NaN, tol=1e-2)
  diff = abs.(vals .- x)
  i = argmin(diff)
  isnan(cell) && return i
  diff[i] <= (0.5 + tol) * abs(cell) ? i : -1 # 在1个网格内, 0.51 cell
end

# cellsize需要是准确的
function findnear((x, y)::Tuple{Real,Real}, lon::AbstractVector, lat::AbstractVector;
  cellx::Real=NaN, celly::Real=NaN, tol=1e-2)

  i = findnear(x, lon; cell=cellx, tol)
  j = findnear(y, lat; cell=celly, tol)
  (i == -1 || j == -1) && (return nothing) # 查找失败
  return i, j
end

findnear(x::Real, y::Real, lon::AbstractVector, lat::AbstractVector; cellx::Real=NaN, celly::Real=NaN, tol=1e-2) =
  findnear((x, y), lon, lat; cellx, celly, tol)

function st_location_exact2(lon::AbstractVector, lat::AbstractVector, points::Vector{Tuple{T,T}};
  rm_empty::Bool=false, cellsize=nothing, tol=1e-2) where {T<:Real}
  isnothing(cellsize) && (cellsize = st_cellsize(lon, lat))
  cellx, celly = cellsize

  locs = map(p -> findnear(p, lon, lat; cellx, celly, tol), points)
  return rm_empty ? _rm_empty(locs) : locs
end

export st_location_exact2
