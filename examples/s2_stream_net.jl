include("main_vis.jl")
using Ipaper, Ipaper.sf, ArchGDAL, DataFrames, RiverGraphs, Test
using Shapefile, DataFrames

pour = Shapefile.Table("data/shp/Pour_十堰_sp8.shp") |> DataFrame
sites = pour.name
points = map(x -> (x.x, x.y), pour.geometry) #|> x -> cat(x..., dims=1)


begin
  g = RiverGraph("./data/十堰_500m_flowdir.tif")
  index_pit = find_pits(g, points)

  info = DataFrame(; name=sites, to=index_pit)
  # level = 3
  min_sto = 5 # maximum(strord) - level

  strord = stream_order(g)
  links = stream_link(g, strord; min_sto)
  add_links!(links, index_pit)
  # links[3575] = maximum(links) + 1 # 竹溪

  ## 手动添加links
  basinId = fillnodata_upstream(g, links, 0)

  river, info_node = fillnodata_upriver(g, links, strord; min_sto, nodata=0)
  net = stream_network(info_node)        # 河网结构
  flow_path(g, info_node, strord; level) # add a depth argument

  for (i, pit) in enumerate(index_pit)
    node = info_node |> d -> d.value_next[d.to.==pit] |> unique
    if !isempty(node)
      r = graph_children(net, only(node))
      println(sites[i], "\t : ", r)
      println()
    end
  end
end

# 孤山   : [18, 19] => 20
# 竹溪   : [15, 21] => 30
# 延坝   : Any[16, [15, 21] => 30] => 22
# 县河   : [11, 12, 13, 14] => 17
# 贾家坊 : [[1, 2, 3, 5] => 7, [4, 6] => 8] => 9
# 房县   : [23, 24, 25, 26, 27] => 28

r = graph_children(net, 9)
unlist(r)
# "松柏（二）", "八亩地": 河道结构较为简单, `min_sto = 5`时无法检测到
# _pours = find_outlet(g.graph, g.toposort, strord; min_sto=2) # dead points

## 子流域存在嵌套关系的如何识别？
index_pit = [indexs_pit[1]]
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
index_pit = [index_pit[1]]
imagesc(_lon, _lat, _data)

begin
  strord_2d, links_2d, basinId_2d =
    Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)
  info_link = getInfo_links(g, links_2d)

  ## [i,j] how to inds
  # _links = fillnodata_upstream(g, links, 0)
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
