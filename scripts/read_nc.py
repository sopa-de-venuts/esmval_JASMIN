from netCDF4 import Dataset
import numpy as np

# nom del fitxer on hi ha les dades de sortida del geogrid.exe
filepath = '/esarchive/exp/cnrm-cm6-1/cmip6-ssp126_i1p1/original_files/cmorfiles/tas_Amon_CNRM-CM6-1_ssp126_r1i1p1f2_gr_201501-210012.nc'
#tas_Amon_CNRM-CM5_rcp85_r10i1p1_200601-205512.nc
#tas_Amon_CNRM-CM5_rcp85_r10i1p1_205601-210012.nc
filepath = '/esarchive/exp/CMIP5/rcp85/cnrm-cm5/original_files/tas_Amon_CNRM-CM5_rcp85_r10i1p1_200601-205512.nc'
nc = Dataset(filepath)
print(nc)
