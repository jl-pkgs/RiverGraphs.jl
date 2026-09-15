using SpatialRasterLite, ArchGDAL
using RiverGraphs, Test
import Graphs


@testset "TauDEM flow direction" begin
  tau = reshape(UInt8.(1:8), 2, 4)
  gis = reshape(UInt8[1, 2, 4, 8, 16, 32, 64, 128], 2, 4)
  @test RiverGraphs.tau2gis(tau) == gis
  @test RiverGraphs.gis2tau(gis) == tau
end


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


@testset "reservoir catchment area" begin
  rg = RiverGraph(UInt8[6 0; 5 0]; lon=[0.0, 1000.0], lat=[0.0, 1000.0], nodata=UInt8(0))
  @test upstream_nodes(rg, 2) == [2, 1]
  @test catchment_area(rg, 2) == 2.0
  @test catchment_stats(rg, 2) == (; outlet=2, upstream_nodes=2, area_km2=2.0)
end


@testset "GuanShan reach table" begin
  rg = RiverGraph(path_flowdir_GuanShan)
  elevation = zeros(rg.ngrid)
  for (rank, node) in enumerate(rg.toposort)
    elevation[node] = rg.ngrid - rank
  end

  result = delineate_reaches(rg; min_sto=4, elevation)
  reaches = result.reaches

  @test nrow(reaches) == 17
  @test count(iszero, reaches.downSegId) == 1
  @test all(>(0), reaches.length)
  @test all(>(0), reaches.slope)
  @test all((reaches.downSegId .== 0) .| (reaches.downSegId .> reaches.segId))
  @test count(>(0), result.river_reach.A) == 125
  @test count(>(0), result.hru_id.A) == rg.ngrid
  @test nrow(result.hrus) == nrow(reaches)
  @test result.hrus.HRUid == result.hrus.hruSegId == reaches.segId

  runoff_grid = zeros(size(rg.index_rev)..., 3)
  for t in axes(runoff_grid, 3), node in 1:rg.ngrid
    i, j = Tuple(rg.index[node])
    runoff_grid[i, j, t] = (node + t) * 1e-9
  end
  inputs = aggregate_hru_runoff(rg, runoff_grid;
    min_sto=4, elevation, cell_area_m2=1.0)

  @test size(inputs.runoff) == (17, 3)
  @test sum(inputs.hrus.area) == rg.ngrid
  for t in axes(runoff_grid, 3)
    grid_volume = sum(runoff_grid[index[1], index[2], t] for index in rg.index)
    @test isapprox(sum(inputs.runoff[:, t] .* inputs.hrus.area), grid_volume)
  end
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
