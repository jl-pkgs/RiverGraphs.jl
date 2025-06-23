using Base.Threads: nthreads
export kinwave_set_subdomains

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
function kinwave_set_subdomains(graph, toposort, index_pit, streamorder, min_sto)
  if nthreads() > 1
    # extract basins (per outlet/pit), assign unique basin id
    n_pits = length(index_pit)
    basin = fill(0, length(toposort))
    basin[index_pit] = [1:n_pits;]
    basin_fill = fillnodata_upstream(graph, toposort, basin, 0)

    # pre-allocate the Vector with indices matching the topological order of the
    # complete domain (upstream neighbors are stored at these indices)
    index_toposort = fill(0, length(toposort))
    for (i, j) in enumerate(toposort)
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
      basin = findall(x -> x == i, basin_fill)
      g, vmap = induced_subgraph(graph, basin)
      toposort_b = topological_sort_by_dfs(g)
      streamorder_subbas = streamorder[vmap]

      subbas = stream_link(g, toposort_b, streamorder_subbas, min_sto)
      subbas_fill = fillnodata_upstream(g, toposort_b, subbas, 0)

      n_subbas = max(length(subbas[subbas.>0]), 1)
      @show "k1", n_subbas

      if n_subbas > 1
        graph_subbas = graph_from_nodes(g, subbas, subbas_fill)
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
      for n in 1:length(v_subbas)
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
          subbas_s = findall(x -> x == s, subbas_fill)
          sg, _ = induced_subgraph(g, subbas_s)
          toposort_sg = topological_sort_by_dfs(sg)
          push!(topo_subbas, basin[subbas_s[toposort_sg]])
          push!(indices_subbas, index_toposort[basin[subbas_s[toposort_sg]]])
        end
      else
        push!(topo_subbas, basin[toposort_b])
        push!(indices_subbas, index_toposort[basin[toposort_b]])
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
