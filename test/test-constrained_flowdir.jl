using RiverGraphs, Test

@testset "D8 flow direction" begin
  dem = reshape([3.0, 2.0, 1.0], 3, 1)
  ldd = d8_flowdir(dem)

  @test vec(ldd) == UInt8[6, 6, 5]
  @test ldd_code(CartesianIndex(1, 1), CartesianIndex(2, 2)) == UInt8(3)
  @test_throws ArgumentError ldd_code(CartesianIndex(1, 1), CartesianIndex(3, 1))
end

@testset "Vector line to 8-connected flow path" begin
  lon = collect(0.0:30.0:120.0)
  lat = collect(120.0:-30.0:0.0)
  line = [(0.0, 120.0), (120.0, 0.0)]

  path = rasterize_flowpath(line, lon, lat)
  @test path == CartesianIndex.(1:5, 1:5)

  for k in 1:length(path)-1
    d = path[k+1] - path[k]
    @test max(abs(d[1]), abs(d[2])) == 1
  end
end

@testset "DEM orientation uses full river profile" begin
  path = CartesianIndex.(1:5, 1:5)
  dem = fill(0.0, 5, 5)
  for (k, I) in enumerate(path)
    dem[I] = k
  end

  @test orient_flowpath(path, dem; direction=:dem) == reverse(path)
  @test orient_flowpath(path, dem; direction=:geometry) == path
end

@testset "Force river paths" begin
  ldd = fill(UInt8(5), 3, 3)
  main = [CartesianIndex(2, 2), CartesianIndex(3, 2)]
  trib1 = [CartesianIndex(1, 1), CartesianIndex(2, 2)]
  trib2 = [CartesianIndex(1, 3), CartesianIndex(2, 2)]

  force_flowpaths!(ldd, [trib1, trib2, main]; outlet=:pit)
  @test ldd[1, 1] == UInt8(3)
  @test ldd[1, 3] == UInt8(9)
  @test ldd[2, 2] == UInt8(6)
  @test ldd[3, 2] == UInt8(5)

  divergent = [
    [CartesianIndex(2, 2), CartesianIndex(3, 2)],
    [CartesianIndex(2, 2), CartesianIndex(2, 3)],
  ]
  @test_throws ArgumentError force_flowpaths!(fill(UInt8(5), 3, 3), divergent)
end

@testset "River-constrained flow direction" begin
  # The DEM makes the centre cell drain south, but trusted hydrography forces
  # the river east through three cells.
  dem = [
    9.0 9.0 9.0
    9.0 5.0 1.0
    9.0 4.0 0.0
  ]
  river = [[
    CartesianIndex(1, 2),
    CartesianIndex(2, 2),
    CartesianIndex(3, 2),
  ]]

  ldd = river_constrained_flowdir(dem, river; outlet=:pit)
  @test ldd[1, 2] == UInt8(6)
  @test ldd[2, 2] == UInt8(6)
  @test ldd[3, 2] == UInt8(5)

  rg = RiverGraph(ldd; nodata=UInt8(0))
  @test rg.ngrid == length(dem)
end
