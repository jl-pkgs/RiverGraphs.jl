import SpatRasters: gdal_nodata, read_gdal

## Note: 
# 1. The dimemsion of `dem`: [lon, lat], which is different from ArcGIS. Hence, DIR
# is different.
# 2. `lat` in the reverse order.
const DIR_GIS = UInt8.([1, 2, 4, 8, 16, 32, 64, 128])
const DIR_TAU = UInt8.(1:8)
const DIR_WFLOW = UInt8.([6, 3, 2, 1, 4, 7, 8, 9])

# DIV = [
#   32    64(N) 128
#   16(W) 0     1(E)
#   8     4(S)  2
# ]
# DIV_WFLOW = [
#   7, 8, 9, 
#   4, 5, 6, 
#   1, 2, 3,
# ]

# "Map from PCRaster LDD value to a CartesianIndex"
# const pcr_dir = [
#   CartesianIndex(-1, -1),  # 1
#   CartesianIndex(0, -1),  # 2
#   CartesianIndex(1, -1),  # 3
#   CartesianIndex(-1, 0),  # 4
#   CartesianIndex(0, 0),  # 5
#   CartesianIndex(1, 0),  # 6
#   CartesianIndex(-1, 1),  # 7
#   CartesianIndex(0, 1),  # 8
#   CartesianIndex(1, 1),  # 9
# ]
"lat reversed `pcr_dir`"
const pcr_dir = [
  CartesianIndex(-1, 1),  # 7
  CartesianIndex(0, 1),  # 8
  CartesianIndex(1, 1),  # 9
  CartesianIndex(-1, 0),  # 4
  CartesianIndex(0, 0),  # 5
  CartesianIndex(1, 0),  # 6
  CartesianIndex(-1, -1),  # 1
  CartesianIndex(0, -1),  # 2
  CartesianIndex(1, -1),  # 3
]


# const NODATA = 0x00
# const DY = [0, 1, 1, 1, 0, -1, -1, -1]
# const DX = [1, 1, 0, -1, -1, -1, 0, 1]
# const DIR = [1, 2, 4, 8, 16, 32, 64, 128]
# const DIR_INV = [16, 32, 64, 128, 1, 2, 4, 8]

# DIR = [4, 2, 1, 128, 64, 32, 16, 8]
# DIR_INV = [64, 32, 16, 8, 4, 2, 1, 128]
# DIR = [3, 2, 1, 8, 7, 6, 5, 4]
# DIR_INV = [7, 6, 5, 4, 3, 2, 1, 8]
# const DX = [0, 1, 1, 1, 0, -1, -1, -1]
# const DY = [1, 1, 0, -1, -1, -1, 0, 1]

function gis2tau(A::AbstractArray)
  R = copy(A)
  for i in 1:8
    replace!(R, DIR_GIS[i] => DIR_TAU[i])
  end
  R
end

function tau2gis(A::AbstractArray)
  R = copy(A)
  for i in 1:8
    replace!(R, DIR_TAU[i] => DIR_GIS[i])
  end
  R
end

function gis2wflow(A::AbstractArray)
  R = copy(A)
  for i in 1:8
    replace!(R, DIR_GIS[i] => DIR_WFLOW[i])
  end
  R
end
