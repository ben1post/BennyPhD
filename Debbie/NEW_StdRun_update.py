import numpy as np
import time
import xarray as xr

from SizeModLogspace_Update import SizeModLogspace

Pnum = 3      # specify the nos of phytoplankton size groups
Ynum = 1       # no. of modeling year
mod_name = 'try'

# Runs - loop
#   Nutrient schemes: Oligotrophic (1.0); Mesotrophic (15.0); Eutrophic (50.0)
#   Mixing regimes: Meromictic (just for observation); polymictic; dimictic; monomictic
N0_list = [1, 15, 50]
mix_list = ['constant', 'medium', 'high']
sol = np.zeros((Ynum*365, Pnum + 4, len(N0_list), len(mix_list)))  # time; no. var; N0 schemes; mixing regimes

# Run solution
starttime = time.ctime()        # save it for documentation purposes
print(mod_name, str(Pnum), "start:", starttime)
start = time.time()             # start the timer

for i in range(len(N0_list)):   # nutrient level
    for j in range(len(mix_list)):  # mixing regimes
        sol[:, :, i, j] = SizeModLogspace(freq=mix_list[j], N0=N0_list[i], numP=Pnum, numYears=Ynum).solution

# Complete
end = time.time()  # stop the timer
endtime = time.ctime()  # save it for documentation purposes
print(mod_name, str(Pnum), "complete:", endtime)
runtime = (end - start) / 3600  # calculate run time (hour)
print(mod_name, str(Pnum), "run time(hour):", runtime)
print(mod_name, str(Pnum), "run time(day):", runtime/24)

# Storing output
out = xr.Dataset(
    {'Nut': (['time', 'n0s', 'MixRgm'], sol[:, 0, :, :]),
     'Zoo': (['time', 'ZooSize', 'n0s', 'MixRgm'], sol[:, 1:3, :, :]),
     'Phy': (['time', 'PhySize', 'n0s', 'MixRgm'], sol[:, 4:, :, :])},
    coords={'time': np.arange(Ynum*365),
            'n0s': N0_list,
            'MixRgm': np.arange(len(mix_list)),
            'ZooSize': [5, 200],
            'PhySize': np.logspace(0, 2, Pnum)},
    attrs={'Title': 'This file contains sol output with different nutrient and mixing regimes',
           'No. phytoplankton size group': Pnum,
           'n0s': 'Nutrient regimes = 1, 15, 50',
           'MixRgm': 'Mixing frequencies = constant, medium (4 times/yr), high (12 times/yr)',
           'LWST units': 'Degree celcius',
           'PAR units': 'Ein m$^{-2}$ d$^{-1}$',
           'MLD units': 'Meters',
           'Variable units': 'µmol NL$^{-1}$',
           'Size units': 'µm ESD',
           'Phyto size range': '1 - 100',
           'Zoo size range': '5 - 200',
           'Start run time': starttime,
           'Complete run time': endtime,
           'Total run time (hours)': runtime})
out.to_netcdf('/people/home/sto/SizeMod/ExpNU_case1.nc')  # export
