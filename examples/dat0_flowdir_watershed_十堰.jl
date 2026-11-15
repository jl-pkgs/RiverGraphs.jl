## 测试`st_shrink`：切取十堰flowdir
using RiverGraphs, SpatialRasterLite, ArchGDAL, Shapefile
using GLMakie, MakieLayers
include("main_vis.jl")

f_pour = "/mnt/z/GitHub/cug-hydro/RiverGraphs.jl/data/shp/Pour_十堰_sp8.shp"
f_flowdir_HuBei = "/mnt/z/GitHub/cug-hydro/RiverGraphs.jl/data/Hubei_500m_flowdir.tif"
f_flowdir = "/mnt/z/GitHub/cug-hydro/RiverGraphs.jl/data/十堰_500m_flowdir.tif"

## clip_region: Hubei -> ShiYan
pour = Shapefile.Table(f_pour)
points = map(x -> (x.x, x.y), pour.geometry) #|> x -> cat(x..., dims=1)

rg = RiverGraph(f_flowdir_HuBei)
ra_basin = st_watershed(rg, points)
ra_basin2 = st_shrink(ra_basin; nodata=0, cellsize_target=0.1)

imagesc(ra_basin)
imagesc(ra_basin2)

## 重新保存ShiYan_flowdir
b = ra_basin2.b
r = read_gdal(f_flowdir, b)
r.A[ra_basin2.A .== 0] .= 0 # 流域外设为 nodata

write_gdal(r, f_flowdir)

## 再起提取流域
rg = RiverGraph(f_flowdir)
ra_basin = st_watershed(rg, points)
imagesc(ra_basin)

# 查找流域范围的出口点
strord = stream_order(rg) # 255 as missing
find_outlet(rg.graph, rg.toposort, strord; min_sto=3)
