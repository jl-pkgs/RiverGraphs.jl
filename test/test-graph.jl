using SpatRasters, ArchGDAL
using RiverGraphs, Test

# flowdir, image(A) should looks normal
@testset "RiverGraph stream_net" begin
  g = RiverGraph(path_flowdir_GuanShan)
  @test Matrix(g)[g.index] == g.data

  ## 4级河流
  level = 2
  @time strord, links, basinId = subbasins(g; level)
  strord_2d, links_2d, basinId_2d =
    Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)
  river, info_node = fillnodata_upriver(g, links, strord; level, nodata=0)
  flow_path(g, info_node, strord; level)

  # index is the index of `topo_subbas`
  subbas_order, indices_subbas, topo_subbas =
    kinwave_set_subdomains(g.graph, g.toposort, [1823], strord; level, parallel=true)

  inds = indices_subbas[1]
  @test g.toposort[inds] == topo_subbas[1]

  @test map(length, info_node.index) == [19, 10, 10, 19, 21,
    8, 25, 25, 8, 14, 4]
  @test size(info_node, 1) == 11

  ## 5级河流
  level = 1
  @time strord, links, basinId = subbasins(g; level)
  strord_2d, links_2d, basinId_2d =
    Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)
  river, info_node = fillnodata_upriver(g, links, strord; level, nodata=0)

  flow_path(g, info_node, strord; level)
  @test map(length, info_node.index) == [19, 19]
end
