using SpatialRasterLite, ArchGDAL
using RiverGraphs, Test
import Graphs


@testset "graph_flow reverse index" begin
  A = UInt8[6 0; 5 0]
  inds, index_rev = active_indices(A, UInt8(0))
  ldd = A[inds]

  graph = graph_flow(ldd, inds, index_rev, pcr_dir)
  graph_compat = graph_flow(ldd, inds, pcr_dir)

  @test Graphs.nv(graph) == 2
  @test Graphs.has_edge(graph, 1, 2)
  @test collect(Graphs.edges(graph)) == collect(Graphs.edges(graph_compat))

  # A downstream index outside the raster must be ignored safely.
  A_boundary = reshape(UInt8[4], 1, 1)
  inds_boundary, index_rev_boundary = active_indices(A_boundary, UInt8(0))
  graph_boundary = graph_flow(A_boundary[inds_boundary], inds_boundary,
    index_rev_boundary, pcr_dir)
  @test Graphs.ne(graph_boundary) == 0
end


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
