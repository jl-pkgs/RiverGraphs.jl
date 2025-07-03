using SpatRasters, ArchGDAL, DataFrames, RiverGraphs, Test
using Shapefile

INDIR = "$(@__DIR__)/../" |> abspath

pours = Shapefile.Table("$INDIR/data/shp/Pour_十堰_sp8.shp") |> DataFrame
rg = RiverGraph("$INDIR//data/十堰_500m_flowdir.tif")

# point2index -> index2point -> 
@testset "point2index" begin
  points = st_points(pours)
  points_next = move2next(rg, points)
  points_bak = index2point(rg, point2index(rg, points_next))
  @test point2index(rg, points_next) == point2index(rg, points_bak)
  @test find_pits(rg, points) == point2index(rg, points_next)
end

@testset "st_stream_network" begin
  r = st_stream_network(rg, pours; min_sto=5)
  display(r)
  @test nrow(r.info_node) == 23
end
