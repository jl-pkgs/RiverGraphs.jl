using SpatialRasterLite, ArchGDAL, Shapefile, RiverGraphs
include("main_vis.jl")

## 
f_pour = "/mnt/z/GitHub/cug-hydro/RiverGraphs.jl/data/shp/Pour_十堰_sp8.shp"
f_flowdir = "/mnt/z/GitHub/cug-hydro/RiverGraphs.jl/data/十堰_500m_flowdir.tif"

pours = Shapefile.Table(f_pour) |> DataFrame
points = st_points(pours)

rg = RiverGraph(f_flowdir)
ra_basin, info_node, net_node = st_stream_network!(rg, pours; min_sto=5);

plot(net_node)
plot(ra_basin)
fig = plot(rg, ra_basin, net_node)
save("Figure1_十堰_river_network.png", fig)

outdir = "./Project_ShiYan/data/basins/"
write_subbasins(rg, info_node, pours; outdir)

## 
# begin
#   rg = RiverGraph("./data/十堰_500m_flowdir.tif")
#   ra_basin, info_node, net_node = st_stream_network!(rg, pours[3:3, :]; min_sto=5)
#   plot(net_node)
#   plot(rg, ra_basin, net_node)
# end

# write_gdal(ra, "Guanshan_subbasins.tif")
# gdal_polygonize("Guanshan_subbasins.tif", "data/shp/Guanshan_subbasins.shp")

## subset basins
## TODO: 为links和points添加编号ID
# [1] 松柏（二）  : 13
# [2] 八亩地      : 16
# [3] 孤山        : [1, 12] => 19
# [4] 竹溪        : [7, 8] => 30
# [5] 延坝        : Any[5, [7, 8] => 30] => 24
# [6] 县河        : [3, 11, 14, 18] => 26
# [7] 贾家坊      : [[20, 21] => 22, [2, 4, 9, 15] => 28] => 29
# [8] 房县        : [6, 10, 17, 23, 25] => 27

# "松柏（二）", "八亩地": 河道结构较为简单, `min_sto = 5`时无法检测到
