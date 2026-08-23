using Test
using RiverGraphs, SpatialRasterLite, ArchGDAL

@testset "RiverGraphs CLI" begin
  cli = joinpath(@__DIR__, "rivergraph")
  flowdir = joinpath(@__DIR__, "..", "data", "GuanShan_flwdir.tif")
  mktempdir() do dir
    arc_out = joinpath(dir, "arcgis")
    tau_out = joinpath(dir, "taudem")
    tau_flowdir = joinpath(dir, "taudem.tif")
    ra = rast(flowdir)
    write_gdal(SpatRaster(RiverGraphs.gis2tau(ra.A), ra), tau_flowdir;
      nodata=ra.nodata[1])

    run(`$cli $flowdir --format arcgis --outputs all --outdir $arc_out --min-sto 6`)
    run(`$cli $tau_flowdir --format taudem --outputs order --outdir $tau_out`)

    stem = "GuanShan_flwdir"
    expected = ["$(stem)_network.csv", "$(stem)_stream_link.tif",
      "$(stem)_stream_order.tif", "$(stem)_subbasins.tif"]
    @test sort(readdir(arc_out)) == expected
    @test readline(joinpath(arc_out, "$(stem)_network.csv")) ==
      "from,to,value,value_next,length,n_node,index"
    @test rast(joinpath(arc_out, "$(stem)_stream_order.tif")).A ==
      rast(joinpath(tau_out, "taudem_stream_order.tif")).A
  end
end
