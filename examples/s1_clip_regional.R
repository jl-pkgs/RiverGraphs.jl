pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, 
  sf, sf2, terra
)

ra = rast("data/Hubei_500m_flowdir.tif")

range = c(109.4, 111.1, 31.7, 33.4)
poly = st_rect(range)

r = crop(ra, poly)
writeRaster(r, "flowdir_r.tif")
