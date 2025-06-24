using Base.Threads: nthreads
export kinwave_set_subdomains

"""
    subbasins_order(g, outlet, max_dist)

Group subbasins (down- to upstream), starting at the `outlet` of the basin, per distance (1
to maximum distance `max_dist`) from the `outlet`, based on the directed acyclic graph `g`
of subbasins. Returns grouped subbasins `order`, ordered from `max_dist` (including
subbasins without upstream neighbor and distance < `max_dist`) to distance 0 (`outlet`).
"""
function subbasins_order(g, outlet, max_dist::Int)
  order = Vector{Vector{Int}}(undef, max_dist + 1)
  order[1] = [outlet]
  for i in 1:max_dist
    v = Vector{Int}()
    for n in order[i]
      ups_nodes = inneighbors(g, n)
      if !isempty(ups_nodes)
        append!(v, ups_nodes)
      end
    end
    order[i+1] = v
  end

  # move subbasins without upstream neighbor (headwater) to index [max_dist+1]
  for i in 1:max_dist
    for s in order[i]
      if isempty(inneighbors(g, s))
        append!(order[max_dist+1], s)
        filter!(e -> e ≠ s, order[i])
      end
    end
  end
  return reverse(order)
end

"""
    graph_from_nodes(graph, subbas, subbas_fill)

Extract directed acyclic graph `g` representing the flow network at subbasin level from
`subbas` containing nodes with a unique id representing subbasin outlets, `subbas_fill` with
subbasin ids for the complete domain, and directed acyclic graph `graph` representing the
flow network for each subbasin cell.
"""
function graph_from_nodes(graph, subbas, subbas_fill)
  n = maximum(subbas)
  g = DiGraph(n)
  for i in 1:n
    idx = findall(x -> x == i, subbas)
    ds_idx = outneighbors(graph, only(idx))
    to_node = subbas_fill[ds_idx]
    if !isempty(to_node)
      add_edge!(g, i, only(to_node))
    end
  end
  return g
end


"""
    function kinwave_set_subdomains(graph, toposort, index_pit, streamorder, min_sto)

Setup subdomains for parallel execution (threading) of the kinematic wave calculation.
Subdomains are subbasins based on a minimum stream order `min_sto` (see also `subbasins(g,
streamorder, toposort, min_sto)`). Subbasins are extracted for each basin outlet in Vector
`index_pit`.

# Arguments
- `graph` directed acyclic graph of the kinematic wave domain
- `toposort` topological order of `graph`
- `index_pit` Vector with basin outlets (pits)
- `streamorder` stream order of the kinematic wave domain
- `min_sto` minimum `streamorder` value

# Output
- `subbas_order` grouped subbasin ids (`Vector{Vector{Int}}`) ordered upstream (first index)
  to downstream (last index)
- `indices_subbas` list of indices per subbasin id stored as `Vector{Vector{Int}}`
- `topo_subbas` topological order per subbasin id stored as `Vector{Vector{Int}}`
"""
function kinwave_set_subdomains(graph, toposort, index_pit, streamorder;
  level::Int=2, min_sto=nothing, parallel::Bool=nthreads() > 1)
  min_sto = get_MinSto(streamorder; level, min_sto)

  if parallel
    # extract basins (per outlet/pit), assign unique basin id
    n_pits = length(index_pit)
    basin = fill(0, length(toposort))
    basin[index_pit] = [1:n_pits;]
    basin_fill = fillnodata_upstream(graph, toposort, basin, 0)

    # pre-allocate the Vector with indices matching the topological order of the
    # complete domain (upstream neighbors are stored at these indices)
    index_toposort = fill(0, length(toposort))
    for (i, j) in enumerate(toposort) # (i, val)
      index_toposort[j] = i
    end

    order_subbas = Vector{Vector{Int}}()
    indices_subbas = Vector{Vector{Int}}()
    topo_subbas = Vector{Vector{Int}}()
    index = Vector{Int}()
    total_subbas = 0

    for i in 1:n_pits
      # extract subbasins per basin, make a graph at the subbasin level, calculate the
      # maximum distance of this graph, and group and order the subbasin ids from
      # upstream to downstream
      index_basin = findall(x -> x == i, basin_fill) # this is index
      g, vmap = induced_subgraph(graph, index_basin)
      toposort_b = topological_sort_by_dfs(g)
      streamorder_subbas = streamorder[vmap]

      links = stream_link(g, toposort_b, streamorder_subbas; level)
      links_fill = fillnodata_upstream(g, toposort_b, links, 0)

      n_subbas = max(length(links[links.>0]), 1)

      if n_subbas > 1
        graph_subbas = graph_from_nodes(g, links, links_fill)
        toposort_subbas = topological_sort_by_dfs(graph_subbas)
        dist = Graphs.Experimental.Traversals.distances(
          Graph(graph_subbas),
          toposort_subbas[end],
        )
        max_dist = maximum([dist; 1])
        v_subbas = subbasins_order(graph_subbas, toposort_subbas[end], max_dist)
      else
        v_subbas = [[1]]
      end
      # subbasins need a unique id (in case of multiple basins/outlets in the
      # kinematic wave domain)
      for n in eachindex(v_subbas)
        v_subbas[n] .= v_subbas[n] .+ total_subbas
      end
      total_subbas += n_subbas
      append!(order_subbas, v_subbas)
      append!(index, 1:length(v_subbas))
      # in case of multiple subbasins calculate topological order per subbasin
      # (subgraph of the corresponding basin graph g), and the indices that match the
      # subbasin topological order
      if n_subbas > 1
        for s in 1:n_subbas
          subbas_s = findall(x -> x == s, links_fill)
          sg, _ = induced_subgraph(g, subbas_s)
          toposort_sg = topological_sort_by_dfs(sg)
          inds = index_basin[subbas_s[toposort_sg]]
          push!(topo_subbas, inds)
          push!(indices_subbas, index_toposort[inds])
        end
      else
        inds = index_basin[toposort_b]
        push!(topo_subbas, inds)
        push!(indices_subbas, index_toposort[inds])
      end
    end
    # reduce the order of subbasin ids by merging groups of subbasins that have the same
    # index (multiple basins/outlets in the kinematic wave domain)
    subbas_order = Vector{Vector{Int}}(undef, maximum(index))
    for m in 1:maximum(index)
      subbas_order[m] = reduce(vcat, order_subbas[index.==m])
    end
  else
    subbas_order = [[1]]
    indices_subbas = [[1:length(toposort);]]
    topo_subbas = [toposort]
  end

  return subbas_order, indices_subbas, topo_subbas
end
