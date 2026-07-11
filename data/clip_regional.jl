using SpatialRasterLite, ArchGDAL

f_in = "/mnt/z/China/ALLChinaRunoff/ChinaBasins/china90_merit/merit90_china_flowdir_arcgis.tif"
f_ref = joinpath(@__DIR__, "Hubei_500m_flowdir.tif")
f_out = joinpath(@__DIR__, "Hubei_90m_flowdir.tif")

b_ref = st_bbox(f_ref)
lon, lat = st_dims(f_in)
ilon = findall(b_ref.xmin .<= lon .<= b_ref.xmax)
ilat = findall(b_ref.ymin .<= lat .<= b_ref.ymax)
b_clip = st_bbox(lon[ilon], lat[ilat])

ra = read_gdal(f_in, b_clip)
write_gdal(ra, f_out)

@info "裁剪完成" f_out size=size(ra.A) bbox=ra.b