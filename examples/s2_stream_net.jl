include("main_vis.jl")
using Ipaper, Ipaper.sf, ArchGDAL, DataFrames, RiverGraphs, Test
using Shapefile, DataFrames

# flowdir, image(A) should looks normal
begin
  g = RiverGraph("./data/十堰_500m_flowdir.tif")
  level = 2
  strord = stream_order(g)

  # min_sto = maximum(strord) - level
  links = stream_link(g, strord; level)
  basinId = fillnodata_upstream(g, links, 0)

  # @time strord, links, basinId = subbasins(g; min_sto)
  strord_2d, links_2d, basinId_2d =
    Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)

  river, info_node = fillnodata_upriver(g, links, 0, strord; level)
  net = stream_network(info_node) # 河网结构
  flow_path(g, info_node, strord; level)
end

d_pour = find_outlet(g.graph, g.toposort, strord; min_sto=2)

pour = Shapefile.Table("data/shp/Pour_十堰_sp8.shp") |> DataFrame
sites = pour.name
points = map(x -> (x.x, x.y), pour.geometry) #|> x -> cat(x..., dims=1)

locs = st_location_exact(g.lon, g.lat, points) # 查找位置
indexs_pit = map(p -> g.index_rev[p[1], p[2]], locs) # 流域出水口的位置
_pours = find_outlet(g.graph, g.toposort, strord; min_sto=2)

## subset_graph
## 子流域存在嵌套关系的如何识别？

function extract_basin(g, index_pit)
  n_pits = length(index_pit)
  basin = fill(0, length(g.toposort))
  
  basin[index_pit] = [1:n_pits;]
  basin_fill = fillnodata_upstream(g.graph, g.toposort, basin, 0)

  basinId_2d = Matrix(g, basin_fill)
  mask = basinId_2d .!== 0
  
  ix, iy = st_shrink(mask, g.lon, g.lat; cellsize_target=0.1)
  _lon, _lat = g.lon[ix], g.lat[iy]
  _data = basinId_2d[ix, iy]
  ## 这是流域形状
  # 重新生成1个graph
end

## 提取流域数据
index_rev = g.index_rev[ix, iy]
index = 




index_pit = [index_pit[1]]
imagesc(_lon, _lat, _data)


## 对于每个流域，提取子流域，并提取河网
# f_flowdir = "data/Hubei_500m_flowdir.tif"
# g = RiverGraph(f_flowdir)
# strord = stream_order(g)
# lon, lat = st_dims(f_flowdir)
## find_root

# n_pits = length(index_pit)
# basin = fill(0, length(g.toposort))
# basin[index_pit] = [1:n_pits;]
# basin_fill = fillnodata_upstream(g.graph, g.toposort, basin, 0)
info_link = getInfo_links(g, links_2d)

begin
  ## [i,j] how to inds
  # xs, ys = get_coord(links_2d .!= 0)
  # vals = filter(x -> x != 0, links_2d)
  _links = fillnodata_upstream(g, links, 0)
  # _links, info_next = link_flow2next(g, links) # 这里是走到了下一个点

  fig = Figure(; size=(1200, 400) .* 1)
  _imagesc!(fig[1, 1], basinId_2d, titles=["BasinId"], fun_axis=rm_ticks!)

  axs, plts, cbar = _imagesc!(fig[1, 2], strord_2d; nodata=-1,
    titles=["Stream Order"], gap=0, fun_axis=rm_ticks!)

  # links, 关键节点
  scatter!(axs[1], info_link.x, info_link.y; color=:black) # colormap not work
  text!(axs[1], info_link.x, info_link.y, text=string.(info_link.link))

  ax = Axis(fig[1, 3])
  plot_graph!(ax, net; offset=0.7)

  colgap!(fig.layout.content[2].content, 10)
  colgap!(fig.layout, 10)
  fig
end
