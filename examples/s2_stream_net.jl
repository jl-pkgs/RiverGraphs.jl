includet("main_vis.jl")
using SpatRasters, ArchGDAL, Shapefile, RiverGraphs, Test
using Printf, RTableTools
import RiverGraphs: isscalar
import Ipaper: table

pours = Shapefile.Table("data/shp/Pour_十堰_sp8.shp") |> DataFrame
points = st_points(pours)

rg = RiverGraph("./data/十堰_500m_flowdir.tif")
ra_basin, info_node, net_node = st_stream_network!(rg, pours; min_sto=5);

plot(net_node)
plot(rg, ra_basin, net_node)

outdir = "./Project_ShiYan/data/basins/"
write_subbasins(rg, info_node, pours; outdir)

# begin
#   rg = RiverGraph("./data/十堰_500m_flowdir.tif")
#   ra_basin, info_node, net_node = st_stream_network!(rg, pours[3:3, :]; min_sto=5)
#   plot(net_node)
#   plot(rg, ra_basin, net_node)
# end

# write_gdal(ra, "Guanshan_subbasins.tif")
# gdal_polygonize("Guanshan_subbasins.tif", "data/shp/Guanshan_subbasins.shp")
plot(net)
plot(_basin)

## subset basins
## TODO: 为links和points添加编号ID
# [1] 松柏（二）  : 1
# [2] 八亩地      : 20
# [3] 孤山        : [9, 10] => 12
# [4] 竹溪        : [8, 13] => 30
# [5] 延坝        : Any[11, [8, 13] => 30] => 14
# [6] 县河        : [15, 16, 17, 18] => 19
# [7] 贾家坊      : [[22, 24] => 25, [21, 23, 26, 27] => 28] => 29
# [8] 房县        : [2, 3, 4, 5, 6] => 7

# "松柏（二）", "八亩地": 河道结构较为简单, `min_sto = 5`时无法检测到
