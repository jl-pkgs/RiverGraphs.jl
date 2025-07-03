using MakieLayers, GLMakie, Graphs, GraphMakie
using Test
using NaNStatistics
using SpatRasters
import MakieLayers: imagesc, imagesc!
import GraphMakie: plot, plot!


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

function plot_discrete!(ax, ra::SpatRaster; nodata=0, kw...)
  lon, lat = st_dims(ra)
  A, colorrange, colors, ticks = build_colorbar(ra.A; nodata)
  imagesc!(ax, lon, lat, A; colors, colorbar=(; ticks), colorrange, force_show_legend=true)
end

function plot_discrete(ra::SpatRaster; nodata=0, kw...)
  lon, lat = st_dims(ra)
  A, colorrange, colors, ticks = build_colorbar(ra.A; nodata)
  imagesc(lon, lat, A; colors, colorbar, kw...)
end


function plot_links!(ax, info_link)
  scatter!(ax, info_link.lon, info_link.lat; color=:black) # colormap not work
  text!(ax, info_link.lon, info_link.lat, text=string.(info_link.link))
end
