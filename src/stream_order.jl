stream_order(g::RiverGraph) = stream_order(g.graph, g.toposort)

"""
    stream_order(g, toposort)

Return the Strahler streamorder based on directed acyclic graph `g` and topological order
`toposort`.
"""
function stream_order(g, toposort)
  n = length(toposort)
  strord = fill(1, n)
  for v in toposort
    inds_up = inneighbors(g, v)
    if !isempty(inds_up)
      sto_up = strord[inds_up]
      if length(findall(x -> x == maximum(sto_up), sto_up)) > 1
        strord[v] = (maximum(sto_up) + 1)
      else
        strord[v] = maximum(sto_up)
      end
    end
  end
  return strord
end
