#!/usr/bin/env julia
# bench_graph.jl — RiverGraph 构建效率测试
#
# 用法: julia --project=. scripts/bench_graph.jl
# 输出: 每阶段耗时(s) + 分配(MiB)，分数据集

using RiverGraphs, SpatialRasterLite, ArchGDAL
import Graphs

const datasets = [
    (; name="Guanshan 2KB", path=joinpath(@__DIR__, "../data/GuanShan_flwdir.tif")),
    (; name="十堰 500m 13KB", path=joinpath(@__DIR__, "../data/十堰_500m_flowdir.tif")),
    (; name="Hubei 500m 1.2MB", path=joinpath(@__DIR__, "../data/Hubei_500m_flowdir.tif")),
    (; name="Hubei 90m 28MB", path=joinpath(@__DIR__, "../data/Hubei_90m_flowdir.tif")),
]

# 带分配统计的计时
function bench(label, f)
    GC.gc()
    before = Base.gc_num()
    t = @elapsed result = f()
    after = Base.gc_num()
    alloc = after.total_allocd - before.total_allocd
    println("  $(rpad(label, 20)) $(rpad(round(t; digits=4), 8))s  $(round(alloc / 1024^2; digits=2))MB")
    result
end

println("RiverGraph 构建基准测试")
println("="^56)

for ds in datasets
    println("\n--- $(ds.name) ---")
    A = read_gdal(ds.path, 1)
    ra = SpatRaster(A, st_bbox(ds.path); nodata=gdal_nodata(ds.path))

    A_gwf = bench("gis2wflow", () -> RiverGraphs.gis2wflow(ra.A))
    inds, rev = bench("active_indices", () -> RiverGraphs.active_indices(A_gwf, ra.nodata[1]))
    ldd = A_gwf[inds]
    graph = bench("graph_flow", () -> graph_flow(ldd, inds, rev, pcr_dir))
    toposort = bench("topo_sort_kahn", () -> topological_sort_kahn(graph))

    nv = Graphs.nv(graph); ne = Graphs.ne(graph)
    println("  顶点: $nv  边: $ne  有环: $(nv - length(toposort) != 0)")
end