library(data.table)
library(ggplot2)
library(sf)

args <- commandArgs(trailingOnly = TRUE)
data_dir <- args[1]
figure_file <- args[2]

reaches <- fread(file.path(data_dir, "reach_table.csv"))
cells <- fread(file.path(data_dir, "hru_cells.csv"))
lines <- fread(file.path(data_dir, "reach_lines.csv"))

parts <- split(lines, lines$segId)
geometries <- lapply(parts, function(x) st_linestring(as.matrix(x[, .(lon, lat)])))
ids <- as.integer(names(parts))
reach_sf <- st_sf(
  reaches[match(ids, segId)],
  geometry = st_sfc(geometries, crs = 4326)
)
st_write(reach_sf, file.path(data_dir, "reaches.gpkg"), delete_dsn = TRUE, quiet = TRUE)

labels <- lines[, .SD[ceiling(.N / 2)], by = segId, .SDcols = c("lon", "lat")]
p <- ggplot(cells, aes(lon, lat)) +
  geom_raster(aes(fill = factor(hruId))) +
  geom_path(
    data = lines,
    aes(group = segId),
    colour = "#132B43",
    linewidth = 0.8
  ) +
  geom_label(
    data = labels,
    aes(label = segId),
    size = 2.6,
    label.size = 0,
    fill = "white"
  ) +
  scale_fill_viridis_d(name = "HRU ID") +
  coord_equal(expand = FALSE) +
  labs(
    title = "GuanShan Reach and HRU Distribution",
    subtitle = sprintf("%d reaches / %d HRUs", nrow(reaches), uniqueN(cells$hruId)),
    x = NULL,
    y = NULL
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(figure_file, p, width = 9, height = 7, dpi = 300, bg = "white")
