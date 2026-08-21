using RiverGraphs, Test

@testset "Hydro-enforced flow direction" begin
  dem = fill(5.0f0, 5, 5)
  river = [[
    CartesianIndex(1, 3),
    CartesianIndex(2, 3),
    CartesianIndex(3, 3),
  ]]

  # Deliberately place a false depression away from the trusted channel.
  sink = CartesianIndex(4, 4)
  dem[sink] = 0.0f0
  dem[CartesianIndex(1, 3)] = 3.0f0
  dem[CartesianIndex(2, 3)] = 2.0f0
  dem[CartesianIndex(3, 3)] = 1.0f0

  state = hydro_enforced_flowdir(dem, river;
    cellsize=(30.0, 30.0),
    bank_drop=0.01,
    min_slope=1e-3,
    boundary_outlets=false,
    outlet=:pit,
    return_dem=true)

  @test state.dem[sink] > dem[sink]
  @test state.flowdir[1, 3] == ldd_code(CartesianIndex(1, 3), CartesianIndex(2, 3))
  @test state.flowdir[2, 3] == ldd_code(CartesianIndex(2, 3), CartesianIndex(3, 3))
  @test state.flowdir[3, 3] == UInt8(5)

  river_cells = Set(state.paths[1])
  for I in CartesianIndices(dem)
    I in river_cells && continue
    @test state.flowdir[I] != UInt8(5)
  end

  rg = RiverGraph(state.flowdir; nodata=UInt8(0))
  @test rg.ngrid == length(dem)
end
