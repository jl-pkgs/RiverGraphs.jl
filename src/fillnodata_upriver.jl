"""
    fillnodata_upriver(rg::RiverGraph, links, streamorder;
        level::Int=2, min_sto=nothing, nodata::Int=0)

与`fillnodata_upbasin`类似，但只填充河道的部分，返回的信息更加详细
"""
function fillnodata_upriver(rg::RiverGraph, links::AbstractVector, strord::AbstractVector;
  level::Int=2, min_sto=nothing, nodata::Int=0)

  g = rg.graph
  toposort = rg.toposort
  min_sto = get_MinSto(strord; level, min_sto)

  ids = filter(x -> x > 0, links)
  nodes = findall(x -> x > 0, links)
  find_node(id) = nodes[findfirst(ids .== id)]

  res = copy(links)
  info_node = []
  # info = Dict{Int,Int}()

  for id_from in reverse(toposort)  # down- to upstream
    strord[id_from] < min_sto && continue
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
        push!(info_node, (; from=find_node(from), to=find_node(to),
          value=from, value_next=to))
      end
    end
  end

  info_node = flow_path(rg, DataFrame(info_node), strord; min_sto)
  res, info_node
end
