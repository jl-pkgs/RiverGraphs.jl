using RiverGraphs, Test
import ArchGDAL

@testset "ArchGDAL river vector adapter" begin
  geojson = """
  {
    "type": "FeatureCollection",
    "features": [
      {
        "type": "Feature",
        "properties": {"id": 1},
        "geometry": {
          "type": "LineString",
          "coordinates": [[0.0, 90.0], [30.0, 60.0], [60.0, 30.0]]
        }
      },
      {
        "type": "Feature",
        "properties": {"id": 2},
        "geometry": {
          "type": "MultiLineString",
          "coordinates": [
            [[60.0, 30.0], [90.0, 0.0]],
            [[30.0, 90.0], [60.0, 60.0]]
          ]
        }
      }
    ]
  }
  """

  mktempdir() do dir
    path = joinpath(dir, "rivers.geojson")
    write(path, geojson)

    lines = read_river_lines(path)
    @test length(lines) == 3
    @test lines[1] == [(0.0, 90.0), (30.0, 60.0), (60.0, 30.0)]
    @test lines[2] == [(60.0, 30.0), (90.0, 0.0)]
    @test lines[3] == [(30.0, 90.0), (60.0, 60.0)]

    lon = collect(0.0:30.0:90.0)
    lat = collect(90.0:-30.0:0.0)
    paths = rasterize_flowpaths(lines, lon, lat)
    @test all(length(path) >= 2 for path in paths)
    @test all(maximum(max(abs(d[1]), abs(d[2]))
      for d in diff(path)) == 1 for path in paths)

    dem = [100.0f0 - i - j for i in eachindex(lon), j in eachindex(lat)]
    state = hydro_enforced_flowdir(dem, path, lon, lat;
      cellsize=(30.0, 30.0), boundary_outlets=true, return_dem=true)
    @test size(state.flowdir) == size(dem)
    @test length(state.paths) == 3
    @test RiverGraph(state.flowdir; nodata=UInt8(0)).ngrid == length(dem)
  end
end
