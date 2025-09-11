using GLMakie, MakieLayers
using SpatialRasterLite, ArchGDAL
using Ipaper

import MakieLayers: imagesc!

function imagesc!(fig, ra::SpatRaster; kw...)
  lon, lat = st_dims(ra)
  imagesc!(fig, lon, lat, ra.A; kw...)
end

function text_center!(ax, ra::SpatRaster)
  A = ra.A[:]
  coords = st_coords(ra)

  ids = setdiff(keys(table(A)), 0) |> collect

  map(id -> begin
      inds = findall(A .== id)
      x = mean(coords[inds, 1])
      y = mean(coords[inds, 2])
      p = (x, y)
      id
      text!(ax, x, y, text="$id"; fontsize=14)
    end, ids)
end


fs = dir("./Project_ShiYan/data/basins/", r"tif$")
sites = str_extract.(basename.(fs), "(?<=_).*(?=_)")

begin
  mar = 0.02
  ncol = 4
  colorrange = (0, 30)
  colors = [nan_color; resample_colors(amwg256, 30)]

  fig = Figure(; size=(1400, 800))
  for k in 1:8
    i = floor(Int, (k - 1) / ncol) + 1
    j = mod(k - 1, ncol) + 1
    
    ra = rast(fs[k])
    ax, plt = imagesc!(fig[i, j], ra; colorrange, colors)
    text!(ax, 0 + mar, 1 - mar, text=sites[k]; font="SimHei", fontsize=16, align=(-0.1, 1.1), space=:relative)
    text_center!(ax, ra)
  end
  fig
  save("Figure1_basinId_河网结构.png", fig)
end
