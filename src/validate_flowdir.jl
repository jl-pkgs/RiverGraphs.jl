"""
    validate_flowdir(ldd; nodata=0)

Validate a D8/PCRaster LDD raster without constructing a `Graphs.jl` graph.

Because D8 is a functional graph (each active cell has at most one downstream
cell), cycle detection only needs a `UInt8` state raster plus one reusable path
stack. This is substantially cheaper than building forward/backward adjacency
lists for very large 30 m domains.

Valid outlets are:
- pit code `5`;
- a direction leaving the raster;
- a direction entering nodata.

Returns `true` when valid. Invalid direction codes and directed cycles raise an
error.
"""
function validate_flowdir(ldd::AbstractMatrix{<:Integer}; nodata=0)
  nodata_ldd = convert(eltype(ldd), nodata)
  state = zeros(UInt8, size(ldd)) # 0=unseen, 1=current path, 2=finished
  stack = Int[]
  sizehint!(stack, 1024)

  LI = LinearIndices(ldd)
  CI = CartesianIndices(ldd)

  @inbounds for start in CartesianIndices(ldd)
    isequal(ldd[start], nodata_ldd) && continue
    state[start] == 2 && continue

    empty!(stack)
    I = start

    while true
      checkbounds(Bool, ldd, I) || break
      isequal(ldd[I], nodata_ldd) && break

      s = state[I]
      s == 2 && break
      if s == 1
        cycle_start = findfirst(==(LI[I]), stack)
        cycle = isnothing(cycle_start) ? copy(stack) : stack[cycle_start:end]
        error("directed cycle detected in flow direction at raster cells $(CI[cycle])")
      end

      state[I] = 1
      push!(stack, LI[I])

      dir = Int(ldd[I])
      1 <= dir <= 9 || throw(ArgumentError(
        "invalid LDD code $dir at $I; expected 1:9 or nodata"))
      dir == 5 && break

      J = I + pcr_dir[dir]
      checkbounds(Bool, ldd, J) || break
      isequal(ldd[J], nodata_ldd) && break
      I = J
    end

    for k in stack
      state[k] = 2
    end
  end

  true
end


export validate_flowdir
