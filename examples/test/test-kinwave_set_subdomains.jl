# index is the index of `topo_subbas`
subbas_order, indices_subbas, topo_subbas = 
  kinwave_set_subdomains(g.graph, g.toposort, [1823], strord, 4)

inds = indices_subbas[1]
@test g.toposort[inds] == topo_subbas[1]

begin
  # topo: 是地形顺序, 
  basinId = fill(0, length(strord))
  for (i, inds) in enumerate(topo_subbas) # why topo_subbas?
    basinId[inds] .= i
  end
  basinId_2d = Matrix(g, basinId)
end



# graph = g.graph
# toposort = g.toposort
# index_pit = [1823]
# n_pits = length(index_pit)
# basin = fill(0, length(toposort))
# basin[index_pit] = [1:n_pits;]
# basin_fill = fillnodata_upbasin(graph, toposort, basin, 0)
