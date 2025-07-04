using SpatRasters, ArchGDAL
using RiverGraphs, Test


# flowdir, image(A) should looks normal
@testset "RiverGraph stream_net" begin
  rg = RiverGraph(path_flowdir_GuanShan)
  @test Matrix(rg)[rg.index] == rg.data

  ## 4级河流
  level = 2
  strord = stream_order(rg)
  links = stream_link(rg, strord; level)
  ra_basin = fillnodata_upbasin(rg, links; nodata=0)
  
  river, info_node = fillnodata_upriver(rg, links, strord; level, nodata=0)

  # index is the index of `topo_subbas`
  subbas_order, indices_subbas, topo_subbas =
    kinwave_set_subdomains(rg.graph, rg.toposort, [835], strord; level, parallel=true)
    # 1823, 835
  inds = indices_subbas[1]
  @test rg.toposort[inds] == topo_subbas[1]

  # @test map(length, info_node.index) == [19, 10, 10, 19, 21,
  #   8, 25, 25, 8, 14, 4]
  # @test size(info_node, 1) == 11

  ## 5级河流
  level = 1
  strord = stream_order(rg)
  links = stream_link(rg, strord; level)
  ra_basin = fillnodata_upbasin(rg, links; nodata=0)
  river, info_node = fillnodata_upriver(rg, links, strord; level, nodata=0)

  flow_path(rg, info_node, strord; level)
  @test map(length, info_node.index) == [19, 19]
end
