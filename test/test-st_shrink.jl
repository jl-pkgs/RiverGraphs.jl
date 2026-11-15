using SpatialRasterLite, ArchGDAL, RiverGraphs, Test
using Shapefile

@testset "st_shrink" begin
  INDIR = "$(@__DIR__)/../" |> abspath

  pours = Shapefile.Table("$INDIR/data/shp/Pour_十堰_sp8.shp") |> DataFrame
  points = st_points(pours) # 范围略小1个网格

  rg = RiverGraph("$INDIR//data/十堰_500m_flowdir.tif")

  ra_basin = st_watershed(rg, points) # 500m, 1/240
  ra_basin2 = st_shrink(ra_basin; nodata=0, cellsize_target=0.05)

  @test bbox2vec(st_bbox(ra_basin2)) ≈ [109.45, 31.7, 111.05, 33.4] atol=1e-10 rtol=0
end
