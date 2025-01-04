#
set echo
#
# rename Band1  to elevation for compatibility with CDF_GEBCO
#
ncrename -v "z,elevation" ETOPO_2022_v1_15s.nc
