"""
    flowdir_qc(ldd, paths; nodata=0)

Summarize the topology quality of a river-constrained D8/LDD raster.

The check rebuilds a `RiverGraph`, so directed cycles are rejected by the
existing topological-sort validation. It then reports:

- exact agreement between channel LDD and the trusted raster river paths;
- total drainage sinks;
- sinks on trusted river cells;
- sinks on the active-domain boundary (outer raster edge or adjacent nodata);
- interior non-river sinks, which normally indicate unresolved depressions or
  broken routing.

For a hydro-enforced basin with a trusted outlet, the main diagnostics should be
`river_match == 1.0` and `interior_nonriver_sinks == 0`.
"""
function flowdir_qc(ldd::AbstractMatrix{<:Integer}, paths::AbstractVector;
  nodata=0)
  nodata_ldd = convert(eltype(ldd), nodata)
  rg = RiverGraph(ldd; nodata=nodata_ldd) # also validates acyclicity
  downstream, river_cells = _river_downstream(paths)

  matched = 0
  for (from, to) in downstream
    checkbounds(Bool, ldd, from) || throw(BoundsError(ldd, from))
    ldd[from] == ldd_code(from, to) && (matched += 1)
  end
  nriver_edges = length(downstream)
  river_match = nriver_edges == 0 ? 1.0 : matched / nriver_edges

  sinks = 0
  river_sinks = 0
  boundary_sinks = 0
  interior_nonriver_sinks = 0
  pits = 0

  @inbounds for v in 1:rg.ngrid
    I = rg.index[v]
    ldd[I] == 5 && (pits += 1)
    isempty(rg.graph.fadjlist[v]) || continue

    sinks += 1
    if I in river_cells
      river_sinks += 1
    elseif _active_domain_boundary(I, ldd, nodata_ldd)
      boundary_sinks += 1
    else
      interior_nonriver_sinks += 1
    end
  end

  (;
    active_cells=rg.ngrid,
    river_cells=length(river_cells),
    river_edges=nriver_edges,
    matched_river_edges=matched,
    river_match,
    sinks,
    river_sinks,
    boundary_sinks,
    interior_nonriver_sinks,
    pits,
  )
end

function _active_domain_boundary(I::CartesianIndex{2},
  ldd::AbstractMatrix, nodata)
  for di in -1:1, dj in -1:1
    (di == 0 && dj == 0) && continue
    J = I + CartesianIndex(di, dj)
    checkbounds(Bool, ldd, J) || return true
    isequal(ldd[J], nodata) && return true
  end
  false
end


export flowdir_qc
