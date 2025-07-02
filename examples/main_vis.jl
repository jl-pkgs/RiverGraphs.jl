using MakieLayers, GLMakie, Graphs, GraphMakie
using Test
import MakieLayers: imagesc
using NaNStatistics
using SpatRasters


function _imagesc!(fig, A; nodata=0, kw...)
  ncol = maximum(A)
  A = replace(A * 1.0, nodata => NaN)
  nlon, nlat = size(A)
  A = reshape(A, nlon, nlat, 1)

  _colors = resample_colors(amwg256, ncol)
  colors = cgrad(_colors, ncol, categorical=true)

  _ticks = 1:ncol
  ticks = _ticks, string.(_ticks)

  imagesc!(fig, 1:nlon, 1:nlat, A; colors,
    colorrange=(0.5, ncol + 0.5),
    force_show_legend=false, colorbar=(; ticks), kw...)
end

function plot_graph!(ax, g; offset=0.2)
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

function plot_graph(g)
  fig = Figure(; size=(800, 600))
  ax = Axis(fig[1, 1])
  plot_graph!(ax, g; offset=0.7)
  fig
end

# plot_graph(net)

function imagesc(ra::SpatRaster)
  lon, lat = st_dims(ra)
  imagesc(lon, lat, ra.A)
end

function imagesc(ax, ra::SpatRaster)
  lon, lat = st_dims(ra)
  imagesc!(ax, lon, lat, ra.A)
end



function build_colorbar(A)
  ncols = nanmaximum(A) |> Int # colors
  _colors = resample_colors(amwg256, ncols)
  colors = cgrad(_colors, ncols, categorical=true)

  _ticks = 1:ncols
  ticks = _ticks, string.(_ticks)
  colorrange = (1, ncols) .+ (-0.5, 0.5)
  colorrange, colors, ticks
end

function plot_basin!(ax, ra::SpatRaster; nodata=0, kw...)
  lon, lat = st_dims(ra)
  A = replace(ra.A, nodata => NaN)
  colorrange, colors, ticks = build_colorbar(A)

  imagesc!(ax, lon, lat, A; colors, colorbar=(;ticks), colorrange, force_show_legend=true)
end

function plot_basin(ra::SpatRaster; nodata=0, kw...)
  lon, lat = st_dims(ra)
  A = replace(ra.A, nodata => NaN)
  colorrange, colors, ticks = build_colorbar(A)

  imagesc(lon, lat, A; colors, colorbar, kw...)
end
