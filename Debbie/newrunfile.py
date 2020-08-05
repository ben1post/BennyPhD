import os
import numpy as np
import time
import xarray as xr
import scipy.interpolate as intrp

from SizeModLogspace_Update import SizeModLogspace

Pnum = 3      # specify the nos of phytoplankton size groups
Ynum = 2       # no. of modeling year

# Import forcing data
os.chdir('/Users/szewing/Desktop/PhD_work/SizeMod/Forcing/Data')
#os.chdir('/people/home/sto/SizeMod/Forcing')

# Temperature and PAR data
LWST = np.genfromtxt('LWST_13p.csv', delimiter=',')
nSSI = np.genfromtxt('nSSI_13p.csv', delimiter=',')
for a in range(8):  # To synchonized the initial and end forcing value
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
    Thimijan, R.W., and R.D. Heins. 1983. Photometric, radiometric and quantum light units of measure:
    A review of procedures for interconversion. HortScience 18(6): 818-822
    """
    extr = 0.45
    conv = 4.57
    sec2day = 86400
    μmol2mol = 1e6
    coef = extr * conv * sec2day / μmol2mol  # 0.1776816
    return netSSI * coef

PAR = nSSI2PAR(nSSI)
PAR[:, 0] = nSSI[:, 0]

# Creating forcing interpolation function

def LWST_int(t):
    spl = intrp.CubicSpline(LWST[:, 0], LWST[:, 5])
    NewT = np.mod(t, 365.)
    lwst = spl(NewT)
    return lwst

def PAR_int(t):
    spl = intrp.CubicSpline(PAR[:, 0], PAR[:, 5])
    NewT = np.mod(t, 365.)
    par = spl(NewT)
    return par

import matplotlib.pyplot as plt
plt.figure()
plt.plot(np.arange(365), LWST_int(np.arange(365)))
plt.plot(np.arange(365), PAR_int(np.arange(365)))

# Sinusoidal mixing forcing
"""
Mix depth of 82.125 is adopted here, the figure is obtained from the mean mixing depth from 8 lakes in Switzerland
     Lake name          Code    Zmix    ThermoZ
# 1. Greifensee         GR      30.0    4.30
# 2. Lake Zürich        LZ      135.0   5.80    
# 3. Hallwilersee       HA      48.0    2.15
# 4. Sempachersee       SE      87.0    0.5
# 5. Baldeggersee       BA      67.0    2.16
# 6. Vierwaldstätersee  VWS/LU  110.0   3.47
# 7. Upper Zürich       UZ      36.0    0.5
# 8. Walensee           WA      144.0   0.5
                                82.125  2.42
"""
# To simply value for modeling: Mixing depth = 82m; thermocline depth = 2.5m

def MLD_Sinu(t, freq):
    """
    https://math.stackexchange.com/questions/100655/cosine-esque-function-with-flat-peaks-and-valleys
    The function produces interpolation for the theoretical sinusoidal mixing depth
    :param t: the time step used by ODEINT solver
    :param freq: the mixing frequency
    :return: an interpolated mixing function
    """
    Zmix = 80       # mixing depth [meter]
    ThermoZ = 2.5   # thermocline depth [meter]
    NewT = np.mod(t, 365.)
    x = np.arange(365)
    if freq == 'constant':
        spl = intrp.CubicSpline(x, np.full(len(x), Zmix))
        mld = spl(NewT)
        return mld
    elif freq == 'medium':
        mld = (-((Zmix - ThermoZ) / 2 + ThermoZ) + (Zmix - ThermoZ) / 2 * np.cos(x / 14.525)) * -1
        spl = intrp.CubicSpline(x, mld)
        mld = spl(NewT)
        return mld
    elif freq == 'high':
        mld = (-((Zmix - ThermoZ) / 2 + ThermoZ) + (Zmix - ThermoZ) / 2 * np.cos(x / 4.825)) * -1
        spl = intrp.CubicSpline(x, mld)
        mld = spl(NewT)
        return mld


def dmdt(t, freq):
    Zmix = 80       # mixing depth [meter]
    ThermoZ = 2.5   # thermocline depth [meter]
    NewT = np.mod(t, 365.)
    x = np.arange(365)
    if freq == 'constant':
        spl = intrp.CubicSpline(x, np.full(len(x), Zmix))
        #mld = spl(NewT)
        dm_spl = spl.derivative()
        dm = dm_spl(NewT)
        return dm
    elif freq == 'medium':
        mld = (-((Zmix - ThermoZ) / 2 + ThermoZ) + (Zmix - ThermoZ) / 2 * np.cos(x / 14.525)) * -1
        spl = intrp.CubicSpline(x, mld)
        dm_spl = spl.derivative()
        dm = dm_spl(NewT)
        return dm
    elif freq == 'high':
        mld = (-((Zmix - ThermoZ) / 2 + ThermoZ) + (Zmix - ThermoZ) / 2 * np.cos(x / 4.825)) * -1
        spl = intrp.CubicSpline(x, mld)
        dm_spl = spl.derivative()
        dm = dm_spl(NewT)
        return dm



"""# Fast plot - LWST, PAR, MLD, dmdt
import matplotlib.pyplot as plt
lab = ['Constant', 'Medium (4 events yr$^{-1}$)', 'High (12 events yr$^{-1}$)']
mix = ['constant', 'medium', 'high']

fig, axs = plt.subplots(2, 2, sharex='all')
plt.tight_layout(pad=2.)

axs[0, 0].plot(np.arange(365), LWST_int(np.arange(365)), 'r', label='LWST')
axs[0, 1].plot(np.arange(365), PAR_int(np.arange(365)), 'orange', label='PAR')
axs[0, 0].set_xlim(0, 365)
#axs[0, 0].set_xlabel('Time [Month]', fontsize=12)
axs[0, 0].set_ylabel('Water Surface Temperature\n[Degree Celcius]', fontsize=12)
axs[0, 1].set_ylabel('PAR [Einstein m$^{-2}$ d$^{-1}$]')
axs[0, 0].set_xticks(np.linspace(0, 365, 12))
axs[1, 0].set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])

for i in range(len(mix)):
    axs[1, 0].plot(np.arange(365), MLD_Sinu(np.arange(365), mix[i]), label=lab[i])
axs[1, 0].set_ylim(85, 0)
axs[1, 0].set_xlim(0, 365)
axs[1, 0].set_xlabel('Time [Month]', fontsize=12)
axs[1, 0].set_ylabel('Mixing depth [meter]', fontsize=12)
axs[1, 0].set_xticks(np.linspace(0, 365, 12))
axs[1, 0].set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
axs[1, 0].legend(loc=1, bbox_to_anchor=(0.85, 1.12), ncol=3)

for i in range(len(mix)):
    axs[1, 1].plot(np.arange(365), dmdt(np.arange(365), mix[i]), label=lab[i])
#axs.set_ylim(85, 0)
axs[1, 1].set_xlim(0, 365)
axs[1, 1].set_xlabel('Time [Month]', fontsize=12)
axs[1, 1].set_ylabel('dMdt [meter]', fontsize=12)
axs[1, 1].set_xticks(np.linspace(0, 365, 12))
axs[1, 1].set_xticklabels(['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'])
axs[1, 1].legend(loc=1, bbox_to_anchor=(0.85, 1.12), ncol=3)
"""

# Runs - loop
#   Nutrient schemes: Oligotrophic (1.0); Mesotrophic (15.0); Eutrophic (50.0)
#   Mixing regimes: Meromictic (just for observation); polymictic; dimictic; monomictic
N0_list = [1, 15, 50]
mix_list = ['constant', 'medium', 'high']
sol = np.zeros((Ynum*365, Pnum + 4, len(N0_list), len(mix_list)))  # time; no. var; N0 schemes; mixing regimes

# Run solution
starttime = time.ctime()        # save it for documentation purposes
print("out", str(Pnum), "start:", starttime)
start = time.time()             # start the timer

for i in range(len(N0_list)):   # nutrient level
    for j in range(len(mix_list)):  # mixing regimes
        sol[:, :, i, j] = SizeModLogspace(mld=MLD_Sinu(t, mix_list[j]), par=PAR_int(t),
                                          sst=LWST_int(t), dmdt=dmdt(t, mix_list[j]),
                                          N0=N0_list[i], numP=Pnum, numYears=Ynum).solution

# Complete
end = time.time()  # stop the timer
endtime = time.ctime()  # save it for documentation purposes
print("out", str(Pnum), "complete:", endtime)
runtime = (end - start) / 3600  # calculate run time (hour)
print("out", str(Pnum), "run time(hour):", runtime)
print("out", str(Pnum), "run time(day):", runtime/24)


# # Storing output
# out = xr.Dataset(
#     {'Nut': (['time', 'n0s', 'MixRgm'], sol[:, 0, :, :]),
#      'Zoo': (['time', 'ZooSize', 'n0s', 'MixRgm'], sol[:, 1:3, :, :]),
#      'Phy': (['time', 'PhySize', 'n0s', 'MixRgm'], sol[:, 4:, :, :])},
#     coords={'time': np.arange(Ynum*365),
#             'n0s': N0_list,
#             'MixRgm': np.arange(numMixRgm),
#             'ZooSize': [5, 200],
#             'PhySize': np.logspace(0, 2, Pnum)},
#     attrs={'Title': 'This file contains sol output with different nutrient and mixing regimes',
#            'No. phytoplankton size group': Pnum,
#            'n0s': 'Nutrient regimes = 1, 15, 50',
#            'LWST units': 'Degree celcius',
#            'PAR units': 'Ein m$^{-2}$ d$^{-1}$',
#            'MLD units': 'Meters',
#            'Variable units': 'µmol NL$^{-1}$',
#            'Size units': 'µm ESD',
#            'Phyto size range': '1 - 100',
#            'Zoo size range': '5 - 200',
#            'Start run time': starttime,
#            'Complete run time': endtime,
#            'Total run time (hours)': runtime})
# out.to_netcdf('/people/home/sto/SizeMod/ExpNU_case1.nc')  # export
