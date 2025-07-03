using MakieLayers, GLMakie, Graphs, GraphMakie
using Test
using NaNStatistics
using SpatRasters, RiverGraphs
import MakieLayers: imagesc, imagesc!
import GraphMakie: plot, plot!

# plot RiverGraph
function plot(rg::RiverGraph, ra_basin, net_node; size=(1600, 600))
  ra_order = SpatRaster(rg, rg.strord)
  info_link = link2point(rg)

  fig = Figure(; size)
  plot!(fig[1, 1], ra_basin; title="BasinId", fun_axis=rm_ticks!)

  ax, plt = plot!(fig[1, 2], ra_order; title="Stream Order", fun_axis=rm_ticks!, nodata=1)
  plot_links!(ax, info_link)

  plot!(Axis(fig[1, 3]), net_node; offset=0.7)

  colgap!(fig.layout.content[2].content, 10)
  colgap!(fig.layout, 10)
  fig
end

# plot graph
function plot!(ax::Axis, g::SimpleDiGraph; offset=0.2)
  p = graphplot!(ax, g;
    nlabels=repr.(1:nv(g)),
    nlabels_color=:blue,
    nlabels_align=(:right, :top),
    arrow_shift=:end, node_size=20, arrow_size=14)
  hidedecorations!(ax)
  hidespines!(ax)

  n = length(p[:node_pos][])
  offsets = map(i -> Point2f(1, 1) .* offset, 1:n)
  p.nlabels_offset[] = offsets
  autolimits!(ax)
end

function plot(g::SimpleDiGraph)
  fig = Figure(; size=(800, 600))
  ax = Axis(fig[1, 1])
  plot!(ax, g; offset=0.7)
  fig
end

# plot
function imagesc(ra::SpatRaster; kw...)
  lon, lat = st_dims(ra)
  imagesc(lon, lat, ra.A; kw...)
end

function imagesc!(ax, ra::SpatRaster; nodata=nothing, kw...)
  A = ra.A[:, :]
  isnothing(nodata) && (A = replace(A, nodata => NaN))
  lon, lat = st_dims(ra)
  imagesc!(ax, lon, lat, A; force_show_legend=false, kw...)
end


## 成熟的函数
function build_colorbar(A; nodata=0)
  A = replace(A, nodata => NaN)

  high = nanmaximum(A) |> Int # colors
  low = nanminimum(A) |> Int
  ncols = high - low + 1

  _colors = resample_colors(amwg256, ncols)
  colors = cgrad(_colors, ncols, categorical=true)

  _ticks = low:high
  ticks = _ticks, string.(_ticks)
  colorrange = (low, high) .+ (-0.5, 0.5)
  A, colorrange, colors, ticks
end


function plot!(ax, ra::SpatRaster{T,N}; nodata=0, kw...) where {T<:Integer,N}
  lon, lat = st_dims(ra)
  A, colorrange, colors, ticks = build_colorbar(ra.A; nodata)
  imagesc!(ax, lon, lat, A; colors, colorbar=(; ticks), colorrange, force_show_legend=true)
end

function plot(ra::SpatRaster{T,N}; nodata=0, kw...) where {T<:Integer,N}
  lon, lat = st_dims(ra)
  A, colorrange, colors, ticks = build_colorbar(ra.A; nodata)
  imagesc(lon, lat, A; colors, colorrange, colorbar=(; ticks), force_show_legend=true, kw...)
end


function plot_links!(ax::Axis, info_link::DataFrame)
  points = info_link.geometry
  link = info_link.link
  scatter!(ax, points; color=:black) # colormap not work
  text!(ax, points, text=string.(link))
end
