using RiverGraphs, SpatRasters, ArchGDAL, Shapefile
using GLMakie, MakieLayers
using SpatRasters: write_gdal
include("main_vis.jl")

pour = Shapefile.Table("data/shp/Pour_十堰_sp8.shp")
points = map(x -> (x.x, x.y), pour.geometry) #|> x -> cat(x..., dims=1)

# f_flowdir = "data/Hubei_500m_flowdir.tif"
f_flowdir = "data/十堰_500m_flowdir.tif"
rg = RiverGraph(f_flowdir)
ra_basin = st_watershed(rg, points)

# strord = stream_order(rg) # 255 as missing
imagesc(ra_basin)

find_outlet(rg.graph, rg.toposort, strord; min_sto=3)
