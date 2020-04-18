import os
import numpy as np
import time
import xarray as xr
from scipy.interpolate import make_interp_spline

from SizeModLogspace import SizeModLogspace

# Import forcing data
os.chdir('/people/home/sto/SizeMod/Forcing')

# Temperature and PAR data
LWST = np.genfromtxt('LWST_13p.csv', delimiter=',')
nSSI = np.genfromtxt('nSSI_13p.csv', delimiter=',')
for a in range(8):      # To synchonized the initial and end forcing value
    LWST[0, a + 1] = LWST[12, a + 1]
    nSSI[0, a + 1] = nSSI[12, a + 1]


# Converting nSSI into PAR
def nSSI2PAR(netSSI):
    """
    This function is to convert nSSI (W m-2) into PAR (Ein m-2 d-1/mol m-2 d-1)
    PAR is a term used to describle radiation in wavelengths usedful for photosynthesis [1]
    PAR is generally accepted to be wavelengths between 400 and 700 nm [1]
    Procedure:
    1. extr = extraction factor (only about 45% of the nSSI is in the 400-700nm range)
    2. conv = W m-2 converts into μmole m-2 s-1 = 4.57 [1,2]
    3. sec2day = converting from s-1 to d-1
    4. μmol2mol = converting from μmol to mol (1 mole = 1 Ein)

    [1] Sager and Farlane, 1997
    Sager, J.C., McFarlane, J.C. 1997. Radiation. In: Plant Growth Chamber Handbook. (Langhans, R.W. and
    Tibbitts, T.W., Eds.) Iowa State University: North Central Regional Research Publication No. 340, Iowa
    Agriculture and Home Economics Experiment Station Special Report no. 99, pp. 1-29.
    [2] Thimijan and Heins, 1983
    Thimijan, R.W., and R.D. Heins. 1983. Photometric, radiometric and quantum light units of measure: A review
    of procedures for interconversion. HortScience 18(6): 818-822
    """
    extr = 0.45
    conv = 4.57
    sec2day = 86400
    μmol2mol = 1e6
    coef = extr * conv * sec2day / μmol2mol  # 0.1776816
    return netSSI * coef


PAR = nSSI2PAR(nSSI)
PAR[:, 0] = nSSI[:, 0]

# Creating forcing curves
LWST_spl = make_interp_spline(LWST[:, 0], LWST[:, 5], k=2)
LWST_spl = LWST_spl(np.arange(365))
PAR_spl = make_interp_spline(PAR[:, 0], PAR[:, 5], k=2)
PAR_spl = PAR_spl(np.arange(365))

# Sinusoidal mixing forcing (Mixing depth = 30m; stratification depth = 1m)
mld_sinu = np.zeros((366, 4))  # days; mixing regimes

mld_sinu[:, 0] = np.full((366,), 30)  # meromictic - no mixing
mld_sinu[:, 1] = (-15.5 + 14.5 * np.cos(np.arange(366) / 9.7)) * -1     # polymictic:hard to have deep mld in polymictic
mld_sinu[:, 2] = (-15.5 + 14.5 * np.cos(np.arange(366) / 29)) * -1      # dimictic
mld_sinu[:, 3] = (-30 + 29 * np.sin(np.arange(366) / 116)) * -1         # monomictic

# define timeline (no. of year to model)
numY = 5
modelT = numY*365

# generate dmdt arrays
dmdt_sinu = np.zeros((365, 4))
for i in range(4):
    dmdt_sinu[:, i] = np.diff(mld_sinu[:, i])  # dmdt (365,)
mld_sinu = mld_sinu[:-1, :]  # (365,)

mld_sinu = np.tile(mld_sinu, (numY, 1))  # multiplicating physical forcing arrays
dmdt_sinu = np.tile(dmdt_sinu, (numY, 1))
LWST_spl = np.tile(LWST_spl, numY)
PAR_spl = np.tile(PAR_spl, numY)

# Runs - loop
# Nutrient schemes: Oligotrophic (1.0); Mesotrophic (15.0); Eutrophic (50.0)
# Mixing regimes: Meromictic (just for observation); polymictic; dimictic; monomictic
sol3 = np.zeros((modelT, 7, 3))  # time; no. var; N0 schemes; mixing regimes
N0_list = [1, 15, 50]

# Run solution
# 3 size classes
starttime = time.ctime()    # save it for documentation purposes
print("out3 start:", starttime)
start = time.time()         # start the timer
for i in range(len(N0_list)):          # nutrient level
    sol3[:, :, i] = SizeModLogspace(mld=mld_sinu[:, 2], par=PAR_spl, sst=LWST_spl, dmdt=dmdt_sinu[:, 2],
                                    N0=N0_list[i], numP=3).solution

# Complete
end = time.time()                       # stop the timer
endtime = time.ctime()                  # save it for documentation purposes
print("out3 complete:", endtime)
runtime = (end - start)/60              # calculate run time
print("out3 run time(min):", runtime)

# Storing output
out3 = xr.Dataset(
    {'Nut': (['time', 'n0s'], sol3[:, 0, :]),
     'Zoo': (['time', 'ZooSize', 'n0s'], sol3[:, 1:3, :]),
     'Phy': (['time', 'PhySize', 'n0s'], sol3[:, 4:, :])},
    coords={'time': np.arange(modelT),
            'n0s': N0_list,
            'ZooSize': [5, 200],
            'PhySize': np.logspace(0, 2, 3)},
    attrs={'Title': 'This file contains sol3 output with different nutrient and mixing regimes',
           'LWST units': 'Degree celcius',
           'PAR units': 'Ein m$^{-2}$ d$^{-1}$',
           'MLD units': 'Meters',
           'Variable units': 'µmol NL$^{-1}$',
           'Size units': 'µm ESD',
           'Phyto size range': '1 - 100',
           'Zoo size range': '5 - 200',
           'Start run time': starttime,
           'Complete run time': endtime,
           'Total run time': runtime})
out3.to_netcdf('/people/home/sto/SizeMod/out3.nc')               # export