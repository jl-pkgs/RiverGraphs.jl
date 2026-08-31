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

@testset "Confluence snapping" begin
  tributary = [CartesianIndex(1, 1), CartesianIndex(2, 1)]
  mainstem = [CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(4, 2)]

  snapped = snap_flowpath_junctions([tributary, mainstem]; maxdist=1)
  @test snapped[1][end] == CartesianIndex(2, 2)
  @test snapped[2] == mainstem

  ldd = fill(UInt8(5), 4, 3)
  force_flowpaths!(ldd, snapped; outlet=:pit)
  @test ldd[2, 1] == ldd_code(CartesianIndex(2, 1), CartesianIndex(2, 2))
  @test ldd[2, 2] == ldd_code(CartesianIndex(2, 2), CartesianIndex(3, 2))
end

@testset "River DEM conditioning" begin
  path = [CartesianIndex(i, 2) for i in 1:5]
  dem = fill(10.0, 5, 3)
  zriver = [12.0, 11.0, 13.0, 10.0, 9.0]
  for (I, z) in zip(path, zriver)
    dem[I] = z
  end

  conditioned = condition_river_dem(dem, [path];
    cellsize=(30.0, 30.0), bank_drop=0.1, min_slope=1e-3)

  # Channel is locally below the 10 m banks and the 13 m DEM bump disappears.
  @test all(conditioned[I] < 10.0 for I in path)
  @test conditioned[path[3]] < dem[path[3]]

  # 0.001 m/m over a cardinal 30 m step gives at least 0.03 m drop.
  @test all(conditioned[path[k]] - conditioned[path[k+1]] >= 0.03 - 1e-10
    for k in 1:length(path)-1)

  # Conditioning is routing-only and never raises terrain.
  @test all(conditioned .<= dem)
end

@testset "Conditioning propagates through confluences" begin
  trib1 = [CartesianIndex(1, 1), CartesianIndex(2, 2)]
  trib2 = [CartesianIndex(1, 3), CartesianIndex(2, 2)]
  mainstem = [CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(4, 2)]

  dem = fill(20.0, 4, 3)
  dem[CartesianIndex(1, 1)] = 15.0
  dem[CartesianIndex(1, 3)] = 14.0
  dem[CartesianIndex(2, 2)] = 16.0 # false high at the confluence
  dem[CartesianIndex(3, 2)] = 13.0
  dem[CartesianIndex(4, 2)] = 12.0

  conditioned = condition_river_dem(dem, [trib1, trib2, mainstem];
    cellsize=(30.0, 30.0), bank_drop=0.0, min_slope=1e-3)

  junction = CartesianIndex(2, 2)
  @test conditioned[junction] < conditioned[CartesianIndex(1, 1)]
  @test conditioned[junction] < conditioned[CartesianIndex(1, 3)]
  @test conditioned[CartesianIndex(3, 2)] < conditioned[junction]
  @test conditioned[CartesianIndex(4, 2)] < conditioned[CartesianIndex(3, 2)]
end

@testset "River-seeded priority flood" begin
  dem = fill(5.0f0, 5, 5)
  seed = CartesianIndex(1, 3)
  dem[seed] = 1.0f0
  sink = CartesianIndex(3, 3)
  dem[sink] = 0.0f0

  flooded = priority_flood_dem(dem, [seed];
    cellsize=(30.0, 30.0), min_slope=1e-3, boundary_outlets=false)

  @test eltype(flooded) == Float32
  @test flooded[seed] == dem[seed]
  @test flooded[sink] > dem[sink]

  ldd = d8_flowdir(flooded; cellsize=(30.0, 30.0))
  for I in CartesianIndices(dem)
    I == seed && continue
    @test ldd[I] != UInt8(5)
  end

  # The overload accepting river paths uses all channel cells as flood seeds.
  river = [[CartesianIndex(1, 3), CartesianIndex(2, 3)]]
  flooded2 = priority_flood_dem(dem, river;
    cellsize=(30.0, 30.0), min_slope=1e-3, boundary_outlets=false)
  @test flooded2[sink] > dem[sink]
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
  @test_throws ArgumentError condition_river_dem(fill(1.0, 3, 3), divergent)
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
