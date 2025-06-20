include("main_vis.jl")
using Ipaper, Ipaper.sf, ArchGDAL, DataFrames, MarrMot

# flowdir, image(A) should looks normal
begin
  f = path_flowdir_GuanShan
  g = RiverGraph(f)

  min_sto = 5
  @time strord, links, basinId = subbasins(g; min_sto)
  strord_2d, links_2d, basinId_2d =
    Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)
  river, info_node = fillnodata_upriver(g, links, 0, strord; min_sto)
  net = stream_network(info_node) # 河网结构
  flow_path(g, info_node, strord; min_sto)
end

@test maximum(info_node.length) ≈ 11.61633647162445

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

begin
  inds = g.index[info_node.index[2]]
  z = deepcopy(links_2d) .* 0
  z[inds] .= 1

  fig = Figure(; size=(600, 400) .* 1)
  _imagesc!(fig[1, 1], links_2d)
  imagesc!(fig[1, 2], z)
  fig
end

begin
  xs, ys = get_coord(links_2d .!= 0)
  vals = filter(x -> x != 0, links_2d)

  # _links = fillnodata_upstream(g, links, 0)
  # _links, info_next = link_flow2next(g, links)
  fig = Figure(; size=(1200, 400).*1)
  _imagesc!(fig[1, 1], basinId_2d, titles=["BasinId"], fun_axis=rm_ticks!)

  axs, plts, cbar = _imagesc!(fig[1, 2], strord_2d; nodata=-1, 
    titles=["Stream Order"], gap=0, fun_axis=rm_ticks!)
  scatter!(axs[1], xs, ys; color=:black) # colormap not work
  text!(axs[1], xs, ys, text=string.(vals))

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
write_gdal(ra, "Guanshan_subbasins.tif")
gdal_polygonize("Guanshan_subbasins.tif", "Guanshan_subbasins.shp")
