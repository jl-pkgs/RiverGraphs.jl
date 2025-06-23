"""
    fillnodata_upstream(g, toposort, data, nodata)

Fill `nodata` upstream cells with the value from the first downstream valid cell, based on
directed acyclic graph `g`, topological order `toposort`, nodata value `nodata` and `data`
containing missing values. Returns filled data `data_out`.
"""
function fillnodata_upstream(g, toposort, data, nodata)
  data_out = copy(data)
  for v in reverse(toposort)  # down- to upstream
    idx_ds = outneighbors(g, v)
    if !isempty(idx_ds)
      if data_out[v] == nodata && data_out[only(idx_ds)] != nodata
        data_out[v] = data_out[only(idx_ds)]
      end
    end
  end
  return data_out
end

fillnodata_upstream(g::RiverGraph, data, nodata) =
  fillnodata_upstream(g.graph, g.toposort, data, nodata)


# 只填充河道的部分
fillnodata_upriver(g::RiverGraph, data, nodata, streamorder; min_sto::Int=4) =
  fillnodata_upriver(g.graph, g.toposort, data, nodata, streamorder; min_sto)

# only for links
function fillnodata_upriver(g::AbstractGraph, toposort,
  links, nodata, streamorder; min_sto::Int=4)

  ids = filter(x -> x > 0, links)
  nodes = findall(x -> x > 0, links)
  find_node(id) = nodes[findfirst(ids .== id)]

  res = copy(links)
  info = []
  # info = Dict{Int,Int}()

  for id_from in reverse(toposort)  # down- to upstream
    streamorder[id_from] < min_sto && continue
    id_to = outneighbors(g, id_from)
    isempty(id_to) && continue

    id_to = only(id_to)
    if res[id_to] != nodata
      if res[id_from] == nodata
        res[id_from] = res[id_to]
      else
        from = res[id_from]
        to = res[id_to]
        # info[from] = to
        push!(info, (; from=find_node(from), to=find_node(to),
          value=from, value_next=to))
      end
    end
  end
  res, DataFrame(info)
end
