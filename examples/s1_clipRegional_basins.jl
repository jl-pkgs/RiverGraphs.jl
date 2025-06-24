using Shapefile, DataFrames
using RiverGraphs
using Ipaper.sf, ArchGDAL
using GLMakie, MakieLayers
using Ipaper: read_flowdir


pour = Shapefile.Table("Z:/GitHub/cug-hydro/Distributed_Hydrology_Forcing/Pour_十堰_sp8.shp")
points = map(x -> (x.x, x.y), pour.geometry) #|> x -> cat(x..., dims=1)

f_flowdir = "data/Hubei_500m_flowdir.tif"
g = RiverGraph(f_flowdir)
strord = stream_order(g)

# lon, lat = st_dims(f_flowdir)
locs = st_location_exact(g.lon, g.lat, points) # 查找位置
index_pit = map(p -> g.index_rev[p[1], p[2]], locs) # 流域出水口的位置

## 1. 提取流域
n_pits = length(index_pit)
basin = fill(0, length(g.toposort))
basin[index_pit] = [1:n_pits;]
basin_fill = fillnodata_upstream(g.graph, g.toposort, basin, 0)

basinId_2d = Matrix(g, basin_fill)
mask = basinId_2d .!== 0
ix, iy = st_shrink(mask, g.lon, g.lat; cellsize_target=0.1)
_lon, _lat = g.lon[ix], g.lat[iy]
_data = basinId_2d[ix, iy]
imagesc(_lon, _lat, _data)


## 提取流域范围
begin
  A = read_gdal(f_flowdir, 1)[:, end:-1:1] # lat was sorted
  A[.!mask] .= UInt8(0) # updateMask

  flowdir = A[ix, iy][:, end:-1:1]
  b = st_bbox(_lon, _lat)
  ra = rast(flowdir, b)
  write_gdal(ra, "data/十堰_500m_flowdir.tif")
end

## 2. 截取目标区域
# lon, lat = st_dims(f)
# lat = reverse(lat)
# gdal_nodata(f_flowdir)
g = RiverGraph("data/十堰_500m_flowdir.tif") # 有失败的
