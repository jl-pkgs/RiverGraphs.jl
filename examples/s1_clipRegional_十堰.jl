using Shapefile, DataFrames
using RiverGraphs
using SpatRasters, ArchGDAL
using SpatRasters: write_gdal
using GLMakie, MakieLayers


pour = Shapefile.Table("data/shp/Pour_十堰_sp8.shp")
points = map(x -> (x.x, x.y), pour.geometry) #|> x -> cat(x..., dims=1)

f_flowdir = "data/Hubei_500m_flowdir.tif"
g = RiverGraph(f_flowdir)

rg = RiverGraph(f_flowdir)
ra_basin = st_watershed(rg, points)
ra_basin2 = st_shrink(basin_fill; nodata=0, cellsize_target=0.1)

imagesc(ra_basin)
imagesc(ra_basin2)
## 
b = ra_basin2.b
r = read_gdal(f_flowdir, b)
fout = "data/十堰_500m_flowdir.tif"
write_gdal(A, fout)

g = RiverGraph(fout) #
