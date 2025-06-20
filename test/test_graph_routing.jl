using Ipaper, Ipaper.sf, ArchGDAL, DataFrames, RiverGraphs, Test


@testset "graph_routing_muskingum!" begin
  # flowdir, image(A) should looks normal
  f = path_flowdir_GuanShan
  g = RiverGraph(f)

  min_sto = 4
  @time strord, links, basinId = subbasins(g; min_sto)
  strord_2d, links_2d, basinId_2d =
    Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)
  river, info_node = fillnodata_upriver(g, links, 0, strord; min_sto)
  flow_path(g, info_node, strord; min_sto)

  net = stream_network(info_node) # 河网结构
  nbasin = nv(net)

  ntime = 100
  data = ones(Float64, ntime, nbasin)
  x = 0.35
  K = 6.0
  l = 5.0
  par_routing = ParamMuskingum(x, K, l)
  graph_routing_muskingum!(net, info_node, data, par_routing)
  @test all(data[:, end] .== 12)
end
