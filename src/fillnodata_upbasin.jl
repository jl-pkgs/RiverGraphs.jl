"""
    fillnodata_upbasin(g, toposort, data, nodata)

Fill `nodata` upstream cells with the value from the first downstream valid
cell, based on directed acyclic graph `g`, topological order `toposort`, nodata
value `nodata` and `data` containing missing values. Returns filled data
`data_out`.
"""
function fillnodata_upbasin(g::AbstractGraph, toposort, data; nodata=0)
  data_out = copy(data)
  for id_from in reverse(toposort)  # down- to upstream
    ids_to = outneighbors(g, id_from) # 
    if !isempty(ids_to)
      if data_out[id_from] == nodata && data_out[only(ids_to)] != nodata
        data_out[id_from] = data_out[only(ids_to)]
      end
    end
  end
  return data_out
end

function fillnodata_upbasin(rg::RiverGraph, data; nodata=0, kw...)
  vec_fill = fillnodata_upbasin(rg.graph, rg.toposort, data; nodata)
  SpatRaster(rg, vec_fill; kw...)
end
