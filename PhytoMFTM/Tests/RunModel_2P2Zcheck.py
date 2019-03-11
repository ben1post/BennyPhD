#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time

# Fitting
from lmfit import minimize, Parameters, Parameter, report_fit

# loading Modules
from PhytoMFTM.ModelClasses import Plankton
import PhytoMFTM.ModelCode as mc


#params.add('xvar', value=0.50, min=0, max=1)
#params.add('yvar', expr='1.0 - xvar')

# set up basic parameters of model:
standardparams = Parameters()
# number of phytoplankton func types
standardparams.add('pfun_num', value=2, vary=False)
# number of zooplankton groups
standardparams.add('zoo_num', value=2, vary=False)
# mld - related
standardparams.add('kappa', value=0.1, vary=False)      # Diffusive mixing across thermocline (m*d^-1)
standardparams.add('deltaD_N', value=0.1, vary=False)   # Nitrate Mineralization rate (d^-1)
standardparams.add('deltaD_Si', value=0.1, vary=False)  # Silicate Mineralization rate (d^-1)

# z - related
standardparams.add('moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
standardparams.add('deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)
#z grazing related
standardparams.add('gr_p', value=0.6, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
standardparams.add('muZ', value=0.1, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# p - related
standardparams.add('kw', value=0.1, vary=False)     # Light attenuation constant (m^-1)
standardparams.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)

standardparams.add('U_N', value=5.0, vary=False)    # Nitrate Half Saturation Constant
standardparams.add('U_Si', value=5.0, vary=False)   # Silicate Half Saturation Constant

standardparams.add('v', value=0.11, vary=False)      # Sinking of Phytoplankton from Mixed Layer

standardparams.add('muP', value=1.6, vary=False)    # Phytoplankton maximum growth rate (d^-1)
standardparams.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)
standardparams.add('Kp', value=0.1, vary=False)     # Zooplankton Grazing assimilation coefficient (-)

standardparams.add('ratioSi', value=1.2, vary=False)  # Silicate ratio


# set up phytoplankton type 1 (e.g. diatoms)
ptype1 = Parameters()
ptype1.add('pt1_U_Si', value=5.0, vary=False)   # Silicate Half Saturation Constant
ptype1.add('pt1_muP', value=1.4, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype1.add('pt1_ratioSi', value=1.25, vary=False)  # Silicate ratio


# set up phytoplankton type 2 (e.g. other)
ptype2 = Parameters()
ptype2.add('pt2_U_Si', value=0., vary=False)     # Silicate Half Saturation Constant
ptype2.add('pt2_muP', value=1.1, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype2.add('pt2_ratioSi', value=0., vary=False)  # Silicate ratio

# set up phytoplankton type 3 (e.g. other)
ptype3 = Parameters()
ptype3.add('pt3_U_Si', value=0., vary=False)     # Silicate Half Saturation Constant
ptype3.add('pt3_U_N', value=1.94, vary=False)    # Nitrate Half Saturation Constant
ptype3.add('pt3_muP', value=0.5, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype3.add('pt3_ratioSi', value=0., vary=False)  # Silicate ratio



# set up zooplankton type 1 (e.g. small zooplankton)
ztype1 = Parameters()
ztype1.add('zt1_gr_p', value=0.6, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
ztype1.add('zt1_muZ', value=0.05, vary=False)    # Zooplankton maximum grazing rate (d^-1)

#ztype1.add('zt1_moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
#ztype1.add('zt1_deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)


# set up zooplankton type 2 (e.g. larger zooplankton)
ztype2 = Parameters()
ztype2.add('zt2_gr_p', value=0.3, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
ztype2.add('zt2_muZ', value=0.1, vary=False)    # Zooplankton maximum grazing rate (d^-1)

#ztype2.add('zt2_moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
#ztype2.add('zt2_deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)







pfn = standardparams['pfun_num'].value
zn = standardparams['zoo_num'].value
print(pfn, zn)

if pfn == 3 and zn == 2:
    print('hi')
    all_params = (standardparams + ptype1 + ptype2 + ptype3 + ztype1 + ztype2)
elif pfn == 2 and zn == 2:
    all_params = (standardparams + ptype1 + ptype2 + ztype1)
    print('hullo')
    all_params = (standardparams + ptype1 + ptype2 + ztype1 + ztype2)
elif pfn == 2 and zn == 1:
    print('cool')
elif pfn == 1 and zn == 1:
    print('oh')
    all_params = (standardparams + ptype1 + ztype1)

# initialize parameters:
N0 = np.mean(mc.NOX)  # Initial Nitrate concentration (mmol*m^-3)
Si0 = np.mean(mc.SiOX)  # Initial Silicate concentration (mmol*m^-3)
Z0 = 0.1 / zn  # Initial Zooplankton concentration (mmol*m^-3)
D0 = 0.01  # Initial Detritus concentration (mmol*m^-3)
P0 = 0.01 / pfn  # Initial Phytoplankton concentration (mmol*m^-3)

initnut = [N0, Si0, D0]
initzoo = [Z0 for i in range(zn)]
initphy = [P0 for i in range(pfn)]
initcond = initnut + initzoo + initphy

timedays_model = np.arange(0., 5 * 365., 1.0)

print(list(all_params)[:])


z = Plankton(all_params,'Zooplankton').init()
p = Plankton(all_params, 'Phytoplankton').init()






#INTEGRATE:
tos = time.time()
print('starting integration')
outarray=odeint(mc.simpleN2P2ZD, initcond, timedays_model, args=(all_params,p,z))
tos1 = time.time()
print('finished after %4.3f sec'%(tos1-tos))





#PLOTTING
a = [plt.plot(outarray[:,3+zn+i]) for i in range(pfn)]
print(pfn)
plt.show()

b = [plt.plot(outarray[:,3+i]) for i in range(zn)]
plt.show()