includet("main_vis.jl")
using SpatRasters, ArchGDAL, DataFrames, RiverGraphs, Test
using Shapefile, DataFrames

pour = Shapefile.Table("data/shp/Pour_十堰_sp8.shp") |> DataFrame
sites = pour.name
points = map(x -> (x.x, x.y), pour.geometry) #|> x -> cat(x..., dims=1)

rg = RiverGraph("./data/十堰_500m_flowdir.tif")
(; lon, lat) = rg
points_next = move2next(rg, _points)

# point2index -> index2point -> 
@testset "point2index" begin
  points_bak = index2point(rg, point2index(rg, points_next))
  @test point2index(rg, points_next) == point2index(rg, points_bak)
  @test find_pits(rg, points) == point2index(rg, points_next)
end


begin
  index_pit = find_pits(rg, points) # 这里是往下移动了一个网格

  info = DataFrame(; name=sites, to=index_pit)
  min_sto = 5 # maximum(strord) - level

  strord = stream_order(rg)
  links = stream_link(rg, strord; min_sto)
  points_link = link2point(rg, links)

  add_links!(links, index_pit)
  # links[3575] = maximum(links) + 1 # 竹溪

  ## 手动添加links
  ra_basin = fillnodata_upstream(rg, links; nodata=0)

  river, info_node = fillnodata_upriver(rg, links, strord; min_sto, nodata=0)
  net = stream_network(info_node)        # 河网结构
  flow_path(rg, info_node, strord; min_sto) # add a depth argument

  for (i, pit) in enumerate(index_pit)
    node = info_node |> d -> d.value_next[d.to.==pit] |> unique
    if !isempty(node)
      r = graph_children(net, only(node))
      println(sites[i], "\t : ", r)
    end
  end
end

plot_graph(net)


# 孤山     : [9, 10] => 12
# 竹溪     : [8, 13] => 30
# 延坝     : Any[11, [8, 13] => 30] => 14
# 县河     : [15, 16, 17, 18] => 19
# 贾家坊   : [[22, 24] => 25, [21, 23, 26, 27] => 28] => 29
# 房县     : [2, 3, 4, 5, 6] => 7
r = graph_children(net, 29)
unlist(r)

# "松柏（二）", "八亩地": 河道结构较为简单, `min_sto = 5`时无法检测到
# _pours = find_outlet(g.graph, g.toposort, strord; min_sto=2) # dead points


begin
  strord_2d, links_2d = Matrix(rg, strord, -1), Matrix(rg, links)
  info_link = getInfo_links(rg, links_2d)

  ## [i,j] how to inds
  # _links = fillnodata_upstream(g, links, 0)
  # _links, info_next = link_flow2next(g, links) # 这里是走到了下一个点

  fig = Figure(; size=(1200, 400) .* 1)
  plot_basin!(fig[1, 1], ra_basin; titles=["BasinId"], fun_axis=rm_ticks!)

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


begin
  fig = Figure(; size=(1200, 800) .* 1)
  plot_basin!(fig[1, 1], ra_basin; titles=["BasinId"], fun_axis=rm_ticks!)
  fig
end


imagesc(ra_basin)

include("main_vis.jl")
plot_basin(ra_basin)
