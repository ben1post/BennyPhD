#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time
import pandas

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
standardparams.add('kappa', value=0.1, min=0.09, max=0.11)      # Diffusive mixing across thermocline (m*d^-1)
standardparams.add('deltaD_N', value=0., vary=False)   # Nitrate Mineralization rate (d^-1)
standardparams.add('deltaD_Si', value=0., vary=False)  # Silicate Mineralization rate (d^-1)

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
standardparams.add('Kp', value=0.5, vary=False)     # Zooplankton Grazing saturation constant (-)

standardparams.add('ratioSi', value=1., vary=False)  # Silicate ratio
"""
add feeding pref params!
Z1P1
Z1P2
Z2P1
Z2P2
"""

# set up phytoplankton type 1 (e.g. diatoms)
#ptype1 = Parameters()
#ptype1.add('pt1_U_Si', value=5.0, vary=False)   # Silicate Half Saturation Constant
#ptype1.add('pt1_muP', value=1.4, vary=False)    # Phytoplankton maximum growth rate (d^-1)
#ptype1.add('pt1_ratioSi', value=1.25, vary=False)  # Silicate ratio


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
#ztype1 = Parameters()
#ztype1.add('zt1_gr_p', value=0.6, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
#ztype1.add('zt1_muZ', value=0.05, vary=False)    # Zooplankton maximum grazing rate (d^-1)

#ztype1.add('zt1_moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
#ztype1.add('zt1_deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)


# set up zooplankton type 2 (e.g. larger zooplankton)
ztype2 = Parameters()
ztype2.add('zt2_gr_p', value=0.6, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
ztype2.add('zt2_muZ', value=0.05, vary=False)    # Zooplankton maximum grazing rate (d^-1)

#ztype2.add('zt2_moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
#ztype2.add('zt2_deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)



all_params = (standardparams)

""" y0 = NRemineralization + NMixing - sum(Gains)  # Nitrate drawdown
    y1 = SiRemineralization + SiMixing - sum(SilicateDrawdown)  # Silicate drawdown
    y2 = sum(ZooMortality) + sum(UnassimilatedProduction) + sum(PhytoMortality) - NRemineralization - SiRemineralization - DetritusMixing   # Detritus

    zoo = [ZooGrowth[i] - ZooMortality[i] - ZooMixing[i] for i in range(zn)]    # Zooplankton losses due to mortality and mixing
    phy = [Gains[i] - PhytoGrazing[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]  # Phytoplankton growth

    outputlist[0] = NRemineralization       - outputlist[0]
    outputlist[1] = NMixing                 - outputlist[1]
    outputlist[2] = SiRemineralization      - outputlist[2]
    outputlist[3] = SiMixing                - outputlist[3]
    outputlist[4] = DetritusMixing          - outputlist[4]

    outputlist[5] = sum(Gains)              - outputlist[5]
    outputlist[6] = sum(SilicateDrawdown)   - outputlist[6]
    outputlist[7] = sum(PhytoMortality)     - outputlist[7]
    outputlist[8] = sum(PhytoGrazing)       - outputlist[8]
    outputlist[9] = sum(PhytoMixing)        - outputlist[9]
    outputlist[10] = sum(PhytoSinking)      - outputlist[10]

    outputlist[11] = sum([LightHarvesting[i] * P[i] for i in range(pfn)])           - outputlist[11]
    outputlist[12] = sum([TemperatureDepGrowth[i] * P[i] for i in range(pfn)])      - outputlist[12]

    outputlist[13] = sum(ZooGrowth)         - outputlist[13]
    outputlist[14] = sum(ZooMortality)      - outputlist[14]
    outputlist[15] = sum(ZooMixing)         - outputlist[15]
    outputlist[16] = sum(UnassimilatedProduction)   - outputlist[16]"""


def setupinitcond(pfn,zn):
    # initialize parameters:
    N0 = 1#np.mean(mc.NOX)  # Initial Nitrate concentration (mmol*m^-3)
    Si0 = 1#np.mean(mc.SiOX)  # Initial Silicate concentration (mmol*m^-3)
    Z0 = 1  # Initial Zooplankton concentration (mmol*m^-3)
    D0 = 1  # Initial Detritus concentration (mmol*m^-3)
    P0 = 1 # Initial Phytoplankton concentration (mmol*m^-3)

    initnut = [N0, Si0, D0]
    initzoo = [Z0 for i in range(zn)]
    initphy = [P0 for i in range(pfn)]
    outputl = [0 for i in range(17)]
    initcond = np.concatenate([initnut, initzoo, initphy])
    #print(type(initcond))
    return initcond


init =setupinitcond(2,2)
time = 1 # days_model = np.arange(0., 5 * 365., 1.0)

z = Plankton(all_params, 'Zooplankton').init()
p = Plankton(all_params, 'Phytoplankton').init()

a1 = mc.phytomftm(init,time,all_params,p,z)
print('2P2Z',a1)

all_params.add('pfun_num', value=2, vary=False)
# number of zooplankton groups
all_params.add('zoo_num', value=1, vary=False)

##########################
init2 =setupinitcond(2,1)

z = Plankton(all_params, 'Zooplankton').init()
p = Plankton(all_params, 'Phytoplankton').init()

a2 = mc.phytomftm(init2,time,all_params,p,z)
print('2P1Z',a2)

all_params.add('pfun_num', value=1, vary=False)
# number of zooplankton groups
all_params.add('zoo_num', value=1, vary=False)
###########################
init3 =setupinitcond(1,1)

z = Plankton(all_params, 'Zooplankton').init()
p = Plankton(all_params, 'Phytoplankton').init()

a3 = mc.phytomftm(init3,time,all_params,p,z)
print('1P1Z',a3)

print('\n a1 phyto', a1[3],a1[4], 'sum phyto', a1[3]+a1[4],
      '\n a1 zoo', a1[5],a1[6], 'sum zoo', a1[5]+a1[6],
      '\n a1 grazing',a1[15],
      '\n a1 zoogrowth', a1[20],
      '\n a1 unassimfeed', a1[23],'\n sum zg uf', a1[20]+a1[23])

print('\n a1 phyto', a1[3],a1[4], 'sum phyto', a1[3]+a1[4],
      '\n a1 zoo', a1[5], 'sum zoo', a1[5],
    '\n a2 grazing',a2[14], '\n a2 zoogrowth', a2[19], '\n a2 unassimfeed', a2[22],'\n  sum zg uf', a2[19]+a2[22])

print('\n a1 phyto', a1[3], 'sum phyto', a1[3],
      '\n a1 zoo', a1[4], 'sum zoo', a1[4],
    '\n a3 grazing',a3[13], '\n a3 zoogrowth', a3[18], '\n a3 unassimfeed', a3[21],'\n  sum zg uf', a3[18]+a3[21])