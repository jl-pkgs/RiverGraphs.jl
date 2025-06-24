include("main_vis.jl")
using Ipaper, Ipaper.sf, ArchGDAL, DataFrames, RiverGraphs, Test

# flowdir, image(A) should looks normal
begin
  f = path_flowdir_GuanShan
  g = RiverGraph(f)

  level = 2
  strord = stream_order(g)
  min_sto = maximum(strord) - level
  links = stream_link(g, strord, min_sto)
  basinId = fillnodata_upstream(g, links, 0)

  # @time strord, links, basinId = subbasins(g; min_sto)
  strord_2d, links_2d, basinId_2d =
    Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)
  river, info_node = fillnodata_upriver(g, links, 0, strord; min_sto)
  net = stream_network(info_node) # 河网结构
  flow_path(g, info_node, strord; min_sto)
end


fig = Figure(; size=(800, 600))
imagesc!(fig, g.lon, g.lat, Matrix(g, g.toposort))
# imagesc!(fig, g.lon, g.lat, basinId_2d)
fig


# @test maximum(info_node.length) ≈ 11.61633647162445

function check_river_length()
  LON, LAT = MarrMot.meshgrid(g.lon, g.lat)
  _inds = info_node.index[1]
  inds = g.index[_inds]
  st_length(LON[inds], LAT[inds])

  fig = Figure(; size=(800, 600))
  ax, plt = imagesc!(fig, g.lon, g.lat, basinId_2d)
  scatter!(ax, LON[inds], LAT[inds])
  fig
end
# check_river_length()

# begin
#   inds = g.index[info_node.index[2]]
#   z = deepcopy(links_2d) .* 0
#   z[inds] .= 1
#   fig = Figure(; size=(600, 400) .* 1)
#   _imagesc!(fig[1, 1], links_2d)
#   imagesc!(fig[1, 2], z)
#   fig
# end


function get_links(g, links_2d)
  con = links_2d .!= 0
  xs, ys = get_coord(con)
  vals = filter(x -> x != 0, links_2d)
  inds = g.index_rev[con][:]
  DataFrame(; x = xs, y = ys, 
    lon = g.lon[xs], lat = g.lat[ys],
    link = vals, index = inds)
end

info_link = get_links(g, links_2d)

begin
  ## [i,j] how to inds
  # xs, ys = get_coord(links_2d .!= 0)
  # vals = filter(x -> x != 0, links_2d)

  # _links = fillnodata_upstream(g, links, 0)
  # _links, info_next = link_flow2next(g, links)
  fig = Figure(; size=(1200, 400).*1)
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


save("Figure1_孤山-河网结构_L5.png", fig; px_per_unit=2)

River = Matrix(g, river)
imagesc(River)

b = st_bbox(path_flowdir_GuanShan)
ra = rast(basinId_2d[:, end:-1:1], b)
# write_gdal(ra, "Guanshan_subbasins.tif")
# gdal_polygonize("Guanshan_subbasins.tif", "data/shp/Guanshan_subbasins.shp")
