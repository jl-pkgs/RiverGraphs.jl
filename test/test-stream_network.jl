using SpatRasters, ArchGDAL, DataFrames, RiverGraphs, Test
using Shapefile

INDIR = "$(@__DIR__)/../" |> abspath

pours = Shapefile.Table("$INDIR/data/shp/Pour_十堰_sp8.shp") |> DataFrame
points = st_points(pours)
rg = RiverGraph("$INDIR//data/十堰_500m_flowdir.tif")

#   # points_next = move2next(rg, points)
#   points_bak = index2point(rg, point2index(rg, points_next))
#   @test point2index(rg, points_next) == point2index(rg, points_bak)
#   @test find_pits(rg, points) == point2index(rg, points_next)
@testset "point2index" begin
  index = point2index(rg, points)
  points_bak = index2point(rg, index)
  @test point2index(rg, points_bak) == index
end

@testset "st_stream_network" begin
  rg = RiverGraph("$INDIR//data/十堰_500m_flowdir.tif")
  ra_basin, info_node, net_node = st_stream_network!(rg, pours; min_sto=5)

  # no warnings
  write_subbasins(rg, info_node, pours[1:1, :]; outdir=".")
  rm("01_松柏（二）_basinId.tif")
  rm("01_松柏（二）_info_node.csv")

  # find_outlet
  out = find_outlet(rg)
  @test out.link == [1, 7, 12, 14, 19, 20, 29]
  @test out.strord == [5, 6, 6, 6, 6, 5, 7]

  display(info_node)
  @test nrow(info_node) == 23
end
