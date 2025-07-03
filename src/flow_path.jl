# 识别每一段的河长
function river_length(river_nodes; lon, lat)
  LON, LAT = meshgrid(lon, lat)
  map(inds -> st_length(LON[inds], LAT[inds]), river_nodes)
end

# 找到两点之间的路径
function flow_path(rg::RiverGraph, from::Int, to::Int, streamorder;
  level::Int=2, min_sto=nothing)
  min_sto = get_MinSto(streamorder; level, min_sto)

  nodes = [from]
  id_from = from

  iter = 0
  while true
    id_to = outneighbors(rg.graph, id_from)
    iter += 1
    if !isempty(id_to)
      id_to = only(id_to)
      id_from = id_to  # 无论如何继续往下走

      streamorder[id_to] < min_sto && continue
      push!(nodes, id_to)
      id_from == to && break
    else
      # 没有下游节点，除非to识别出错，否则不会到达这里
      @error "No downstream node found for $id_from"
      break
    end
  end
  nodes
end

function flow_path(rg::RiverGraph, info_node::DataFrame, streamorder;
  level::Int=2, min_sto=nothing)
  min_sto = get_MinSto(streamorder; level, min_sto)

  index_g = []
  for i in 1:nrow(info_node)
    from, to = info_node[i, [:from, :to]]
    inds = flow_path(rg, from, to, streamorder; level, min_sto)
    push!(index_g, inds)
  end
  n_node = length.(index_g)
  # 河道上的点, 根据subbasins可以获取流域的点
  index = map(I -> rg.index[I], index_g)
  len = river_length(index; lon=rg.lon, lat=rg.lat)
  cbind(info_node; length=len, n_node, index=index_g) # add `length`, `n_node`, `index`
end

"""
links to main streams

# Return
- `info`: 流动方向的记录
"""
function link_flow2next(g::RiverGraph, links::AbstractVector)
  (; graph, toposort) = g
  n = length(toposort)
  links2 = fill(0, n)

  # info = Dict{Int,Int}()
  info = []
  for v in toposort
    links[v] == 0 && continue

    ds_nodes = outneighbors(graph, v)
    next = !isempty(ds_nodes) ? only(ds_nodes) : v

    links2[next] = links[v]
    # info[v] = next
    push!(info, (; from=v, to=next, value=links[v]))
  end
  links2, DataFrame(info)
end

link_flow2next(g::RiverGraph, links_2d::AbstractMatrix) =
  link_flow2next(g, links_2d[g.index])


export link_flow2next
export flow_path
