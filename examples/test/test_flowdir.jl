#!/usr/bin/env julia
# 流向栅格质量检查
# 用法: julia --project=. scripts/inspect_flowdir.jl [path_to_flowdir.tif]

using SpatialRasterLite, ArchGDAL, RiverGraphs
import Graphs

const f_default = joinpath(@__DIR__, "../data/Hubei_500m_flowdir.tif")
const f = isempty(ARGS) ? f_default : ARGS[1]

# 1. 元数据
ds = ArchGDAL.read(f)
band = ArchGDAL.getband(ds, 1)
A = ArchGDAL.read(ds, 1)
n_total = length(A)
nodata = UInt8(ArchGDAL.getnodatavalue(band))
n_nodata = count(==(nodata), A)
n_active = n_total - n_nodata

println("文件: $f")
println("尺寸: $(ArchGDAL.width(ds)) × $(ArchGDAL.height(ds))  ($(eltype(A)))  nodata=$nodata")
println("有效: $n_active / $n_total ($(round(Int, n_active/n_total*100))%)")

# 2. 流向值分布（ArcGIS D8）
println("\n流向值分布:")
for v in sort(unique(A))
    pct = round(Int, count(==(v), A) / n_total * 100)
    println("  $v  ($pct%)")
end

# 3. 构建 RiverGraph
@time rg = RiverGraph(f)

n_pits = count(==(UInt8(5)), rg.data)
outdeg = Graphs.outdegree(rg.graph, 1:rg.ngrid)
n_outlets = count(==(0), outdeg)
println("顶点: $(rg.ngrid)  边: $(Graphs.ne(rg.graph))  pit: $n_pits  出口: $n_outlets")
println("toposort 长度: $(length(rg.toposort))")
# cycle 检测已在 topological_sort_kahn 内部完成，无需二次扫描
