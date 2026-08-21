using RiverGraphs, Test
import SpatialRasterLite, ArchGDAL


include("test-st_shrink.jl")
include("test-graph.jl")
include("test-stream_network.jl")
include("test-constrained_flowdir.jl")
include("test-hydro_enforced_flowdir.jl")
include("test-archgdal-ext.jl")
