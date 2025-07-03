includet("main_vis.jl")
using SpatRasters, ArchGDAL, DataFrames, RiverGraphs, Test
using Shapefile, DataFrames

pours = Shapefile.Table("data/shp/Pour_十堰_sp8.shp") |> DataFrame
sites = pours.name
points = st_points(pours)

rg = RiverGraph("./data/十堰_500m_flowdir.tif")
(; lon, lat) = rg
points_next = move2next(rg, points)

# point2index -> index2point -> 
@testset "point2index" begin
  points_bak = index2point(rg, point2index(rg, points_next))
  @test point2index(rg, points_next) == point2index(rg, points_bak)
  @test find_pits(rg, points) == point2index(rg, points_next)
end

r = st_stream_network(rg, pours; min_sto=5)
plot(r.net_node)
# plot(r.net)

## TODO: 为links和points添加编号ID
# 孤山     : [9, 10] => 12
# 竹溪     : [8, 13] => 30
# 延坝     : Any[11, [8, 13] => 30] => 14
# 县河     : [15, 16, 17, 18] => 19
# 贾家坊   : [[22, 24] => 25, [21, 23, 26, 27] => 28] => 29
# 房县     : [2, 3, 4, 5, 6] => 7
graph_children(r.net_node, 29)
unlist(r)

# "松柏（二）", "八亩地": 河道结构较为简单, `min_sto = 5`时无法检测到
# _pours = find_outlet(g.graph, g.toposort, strord; min_sto=2) # dead points
include("main_vis.jl")

begin
  ra_order = SpatRaster(rg, r.strord)
  ra_link = SpatRaster(rg, r.links)
  info_link = getInfo_links(rg, ra_link.A)

  fig = Figure(; size=(1600, 600) .* 1)
  plot_discrete!(fig[1, 1], ra_basin; title="BasinId", fun_axis=rm_ticks!)

  ax, plt = plot_discrete!(fig[1, 2], ra_order; title="Stream Order", fun_axis=rm_ticks!, nodata=1)
  plot_links!(ax, info_link)

  plot!(Axis(fig[1, 3]), r.net_node; offset=0.7)

  colgap!(fig.layout.content[2].content, 10)
  colgap!(fig.layout, 10)
  fig
end
