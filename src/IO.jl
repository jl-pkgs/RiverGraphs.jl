import Ipaper
import Ipaper.sf: gdal_nodata, read_gdal

## Note: 
# 1. The dimemsion of `dem`: [lon, lat], which is different from ArcGIS. Hence, DIR
# is different.
# 2. `lat` in the reverse order.
const DIR_GIS = UInt8.([1, 2, 4, 8, 16, 32, 64, 128])
const DIR_TAU = UInt8.(1:8)
const DIR_WFLOW = UInt8.([6, 3, 2, 1, 4, 7, 8, 9])
# # 按照这种方法
# const pcr_dir = [
#     CartesianIndex(-1, -1),  # 1, 8
#     CartesianIndex(0, -1),   # 2, 4
#     CartesianIndex(1, -1),   # 3, 2
#     CartesianIndex(-1, 0),   # 4, 16
#     CartesianIndex(0, 0),    # 5, 0
#     CartesianIndex(1, 0),    # 6, 1
#     CartesianIndex(-1, 1),   # 7, 32
#     CartesianIndex(0, 1),    # 8, 64
#     CartesianIndex(1, 1),    # 9, 128
# ]

const NODATA = 0x00
const DY = [0, 1, 1, 1, 0, -1, -1, -1]
const DX = [1, 1, 0, -1, -1, -1, 0, 1]
const DIR = [1, 2, 4, 8, 16, 32, 64, 128]
const DIR_INV = [16, 32, 64, 128, 1, 2, 4, 8]

# DIV = [
#   32 64 128
#   16 0  1
#   8  4  2
# ]

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


# 0~9
function gis2wflow(A::AbstractArray; nodata::UInt8)
  R = copy(A)
  for i in 1:8
    replace!(R, DIR_GIS[i] => DIR_WFLOW[i])
  end
  replace!(R, UInt8(nodata) => UInt8(5))
  R
end

"""
- `pit`: flowdir is `0`
"""
function read_flowdir(f::String)
  A_gis = read_gdal(f)[:, end:-1:1] # 修正颠倒的lat
  nodata = gdal_nodata(f)[1]
  A = gis2wflow(A_gis; nodata) # nodata as pit
  # replace!(A, nodata => 0) # replace missing value with 0
  A
end

export read_flowdir
