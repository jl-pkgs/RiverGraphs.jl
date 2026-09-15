using ArchGDAL, RiverGraphs
using DataFrames: DataFrame, Not, eachrow
using RTableTools: fwrite
using SpatialRasterLite: write_gdal

rg = RiverGraph(path_flowdir_GuanShan)
ntime = 24
runoff_grid = zeros(size(rg.index_rev)..., ntime)
for t in 1:ntime, node in 1:rg.ngrid
  i, j = Tuple(rg.index[node])
  storm = max(0.0, sinpi(t / ntime))
  runoff_grid[i, j, t] = 1e-6 * storm * (0.75 + 0.5 * node / rg.ngrid)
end
result = aggregate_hru_runoff(rg, runoff_grid; min_sto=4)

@assert nrow(result.reaches) == 17
@assert nrow(result.hrus) == 17
@assert count(iszero, result.reaches.downSegId) == 1

outdir = joinpath(@__DIR__, "data", "s3_reach_table")
figfile = joinpath(@__DIR__, "..", "images", "s3_reach_hru.png")
mkpath(outdir)
mkpath(dirname(figfile))

hru = result.hru_id.A[rg.index]
reach = result.river_reach.A[rg.index]
cells = DataFrame(
  lon=[rg.lon[index[1]] for index in rg.index],
  lat=[rg.lat[index[2]] for index in rg.index],
  hruId=hru,
  segId=reach,
)

lines = DataFrame(segId=Int[], vertex=Int[], lon=Float64[], lat=Float64[])
for row in eachrow(result.reaches), (vertex, node) in enumerate(row.index)
  index = rg.index[node]
  push!(lines, (row.segId, vertex, rg.lon[index[1]], rg.lat[index[2]]))
end

runoff = DataFrame(HRUid=result.hrus.HRUid)
for t in axes(result.runoff, 2)
  runoff[!, Symbol("t", lpad(string(t), 3, '0'))] = result.runoff[:, t]
end

fwrite(result.reaches[:, Not(:index)], joinpath(outdir, "reach_table.csv"))
fwrite(result.hrus, joinpath(outdir, "hru_table.csv"))
fwrite(runoff, joinpath(outdir, "runoff_hru.csv"))
fwrite(cells, joinpath(outdir, "hru_cells.csv"))
fwrite(lines, joinpath(outdir, "reach_lines.csv"))
write_gdal(result.river_reach, joinpath(outdir, "river_reach.tif"); nodata=0)
write_gdal(result.hru_id, joinpath(outdir, "hru_id.tif"); nodata=0)

rscript = get(ENV, "RSCRIPT", "Rscript")
run(`$rscript $(joinpath(@__DIR__, "plot_s3_reach_table.R")) $outdir $figfile`)

println("reaches/HRUs: ", nrow(result.reaches), "/", nrow(result.hrus))
println("runoff: ", size(runoff_grid), " grid → ", size(result.runoff), " HRU × time")
println("data: ", outdir)
println("figure: ", figfile)
