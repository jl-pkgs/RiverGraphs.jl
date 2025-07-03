using SpatRasters, ArchGDAL, RiverGraphs, Test
using Shapefile

@testset "st_shrink" begin
  INDIR = "$(@__DIR__)/../" |> abspath

  pours = Shapefile.Table("$INDIR/data/shp/Pour_十堰_sp8.shp") |> DataFrame
  points = st_points(pours) # 范围略小1个网格

  rg = RiverGraph("$INDIR//data/十堰_500m_flowdir.tif")

  ra_basin = st_watershed(rg, points) # 500m, 1/240
  ra_basin2 = st_shrink(ra_basin; nodata=0, cellsize_target=0.05)

  @test st_bbox(ra_basin2) == bbox(109.44999999999996, 31.7000000000001, 111.0499999999985, 33.400000000000006)
end
