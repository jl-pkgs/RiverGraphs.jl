using RiverGraphs, Test

@testset "Flow direction validation" begin
  valid = reshape(UInt8[6, 6, 5], 3, 1)
  @test validate_flowdir(valid)

  cycle = reshape(UInt8[6, 4], 2, 1)
  @test_throws ErrorException validate_flowdir(cycle)

  invalid_code = reshape(UInt8[10], 1, 1)
  @test_throws ArgumentError validate_flowdir(invalid_code)

  # Flow may legally leave the active domain or enter nodata.
  outward = reshape(UInt8[4], 1, 1)
  @test validate_flowdir(outward)

  masked = UInt8[6 0; 5 0]
  @test validate_flowdir(masked)
end
