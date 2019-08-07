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
from PhytoMFTM.ModelClasses import Plankton, Forcing
import PhytoMFTM.ModelCode as mc


#params.add('xvar', value=0.50, min=0, max=1)
#params.add('yvar', expr='1.0 - xvar')

# set up basic parameters of model:
standardparams = Parameters()

# number of phytoplankton func types
# standardparams.add('pfun_num', value=4, vary=False)
# number of zooplankton groups
# standardparams.add('zoo_num', value=2, vary=False)

# mld - related
standardparams.add('kappa', value=0.1, vary=False) # min=0.09, max=0.11) # Diffusive mixing across thermocline (m*d^-1)
standardparams.add('deltaD_N', value=0.05, vary=False)   # Nitrate Remineralization rate (d^-1)

standardparams.add('kw', value=0.04, vary=False)     # Light attenuation constant of water (m^-1)

standardparams.add('kc', value=0.03, vary=False)      # Light attenuation via phytoplankton pigment (m^-1)
standardparams.add('alpha', value=0.15, vary=False)  # initial slope of the P-I curve
standardparams.add('VpMax', value=1., vary=False)    # maximum photosynthetic rate

standardparams.add('v', value=0.01, vary=False)      # Sinking of Phytoplankton from Mixed Layer
standardparams.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)

# p - related
standardparams.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)

standardparams.add('ratioSi', value=0, vary=False)  # Silicate ratio ## FALLBACK PARAM FOR OTHER PFTs
standardparams.add('U_N', value=0, vary=False)    # Nitrate Half Saturation Constant
standardparams.add('U_Si', value=0, vary=False)   # Silicate Half Saturation Constant
standardparams.add('muP', value=0, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 1 (e.g. DIATOMS)
ptype1 = Parameters()
ptype1.add('pt1_ratioSi', value=1, vary=False)  # Silicate ratio
ptype1.add('pt1_U_Si', value=1.5, vary=False)   # Silicate Half Saturation Constant
ptype1.add('pt1_U_N', value=1.5, vary=False)    # Nitrate Half Saturation Constant
ptype1.add('pt1_muP', value=1.2, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 2 (e.g. Haptos)
ptype2 = Parameters()
#ptype2.add('pt2_OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)
ptype2.add('pt2_U_N', value=1., vary=False)    # Nitrate Half Saturation Constant
ptype2.add('pt2_muP', value=1., vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 3 (e.g. Cyanos)
ptype3 = Parameters()
ptype3.add('pt3_U_N', value=0.7, vary=False)    # Nitrate Half Saturation Constant
ptype3.add('pt3_muP', value=0.5, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 4 (e.g. Dinos)
ptype4 = Parameters()
ptype4.add('pt4_U_N', value=0.8, vary=False)    # Nitrate Half Saturation Constant
ptype4.add('pt4_muP', value=0.5, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 5 (e.g. Others)
ptype5 = Parameters()
ptype5.add('pt5_U_N', value=0.8, vary=False)    # Nitrate Half Saturation Constant
ptype5.add('pt5_muP', value=0.8, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# z - related
#z grazing related
standardparams.add('moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
standardparams.add('deltaZ', value=0.75, vary=False)    # Zooplankton Grazing assimilation coefficient (-)
standardparams.add('deltaLambda', value=0.75, vary=False)    # Zooplankton Inter-Grazing assimilation coefficient (-)
standardparams.add('muIntGraze', value=0.001, vary=False)  # InterZooGrazing maximum grazing rate
standardparams.add('kIntGraze', value=1., vary=False)  # InterZooGrazing saturation constant

standardparams.add('Kp', value=0, vary=False)     # Zooplankton Grazing saturation constant (-)
standardparams.add('pred', value=0, vary=False)  # quadratic higher order predation rate on zooplankton
standardparams.add('muZ', value=0, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# set up zooplankton type 1 (e.g. MIKRO zooplankton)
ztype1 = Parameters()
ztype1.add('zt1_muZ', value=0.3, vary=False)    # Zooplankton maximum grazing rate (d^-1)

ztype1.add('zt1_Kp', value=.5, vary=False)       # Zooplankton Grazing saturation constant (-)
ztype1.add('zt1_pred', value=0.01, vary=False)    # quadratic higher order predation rate on zooplankton

# set up zooplankton type 2 (e.g. MESO zooplankton)
ztype2 = Parameters()
ztype2.add('zt2_muZ', value=0.1, vary=False)    # Zooplankton maximum grazing rate (d^-1)

ztype2.add('zt2_Kp', value=.8, vary=False)       # Zooplankton Grazing saturation constant (-)
ztype2.add('zt2_pred', value=0.01, vary=False)    # quadratic higher order predation rate on zooplankton

"""
add feeding pref params!
Z1P1
Z1P2
Z2P1
Z2P2
"""
# MIKRO
ztype1.add('zt1_P1', value=0, vary=False)  # Diatoms
ztype1.add('zt1_P2', value=0.25, vary=False)  # Hapto
ztype1.add('zt1_P3', value=0.25, vary=False)  # Cyano
ztype1.add('zt1_P4', value=0.25, vary=False)  # Dino
ztype1.add('zt1_P5', value=0.25, vary=False)  # Others
# MESO
ztype2.add('zt2_P1', value=0.33, vary=False)
ztype2.add('zt2_P2', value=0.25, vary=False)
ztype2.add('zt2_P3', value=0.125, vary=False)
ztype2.add('zt2_P4', value=0.125, vary=False)
ztype2.add('zt2_P5', value=0.125, vary=False)

# inter zoo feeding
# MIKRO
ztype1.add('zt1_Zint_feed1', value=0, vary=False)
ztype1.add('zt1_Zint_feed2', value=0, vary=False)
# MESO
ztype2.add('zt2_Zint_feed1', value=0.01, vary=False)
ztype2.add('zt2_Zint_feed2', value=0, vary=False)

# CONVERT FEEDPREFS TO GRAZEPREF FOR CALCULATION OF GRAZING
ztype1.add('zt1_Zint_grazed1', value=ztype1['zt1_Zint_feed1'].value, vary=False)
ztype1.add('zt1_Zint_grazed2', value=ztype2['zt2_Zint_feed1'].value, vary=False)

ztype2.add('zt2_Zint_grazed1', value=ztype1['zt1_Zint_feed2'].value, vary=False)
ztype2.add('zt2_Zint_grazed2', value=ztype2['zt2_Zint_feed2'].value, vary=False)

ptype1.add('pt1_Z1', value=ztype1['zt1_P1'].value, vary=False)
ptype1.add('pt1_Z2', value=ztype2['zt2_P1'].value, vary=False)

ptype2.add('pt2_Z1', value=ztype1['zt1_P2'].value, vary=False)
ptype2.add('pt2_Z2', value=ztype2['zt2_P2'].value, vary=False)

ptype3.add('pt3_Z1', value=ztype1['zt1_P3'].value, vary=False)
ptype3.add('pt3_Z2', value=ztype2['zt2_P3'].value, vary=False)

ptype4.add('pt4_Z1', value=ztype1['zt1_P4'].value, vary=False)
ptype4.add('pt4_Z2', value=ztype2['zt2_P4'].value, vary=False)

ptype5.add('pt5_Z1', value=ztype1['zt1_P5'].value, vary=False)
ptype5.add('pt5_Z2', value=ztype2['zt2_P5'].value, vary=False)


def setupinitcond(pfn,zn):
    # initialize parameters:
    N0 = 2  # Initial Nitrate concentration (mmol*m^-3)
    Si0 = 2  # Initial Silicate concentration (mmol*m^-3)
    Z0 = 0.01 / zn  # Initial Zooplankton concentration (mmol*m^-3)
    D0 = 0.0  # Initial Detritus concentration (mmol*m^-3)
    P0 = 0.01 / pfn  # Initial Phytoplankton concentration (mmol*m^-3)

    initnut = [N0, Si0, D0]
    initzoo = [Z0 for i in range(zn)]
    initphy = [P0 for i in range(pfn)]
    outputl = [0 for i in range(20)]
    initcond = np.concatenate([initnut, initzoo, initphy, outputl])
    #print(type(initcond))
    return initcond


def runmodel(all_params, initcond, forcing):
    print(list(all_params)[:])

    z = Plankton(all_params, 'Zooplankton').init()
    p = Plankton(all_params, 'Phytoplankton').init()
    fx = Forcing(forcing)
    # INTEGRATE:
    tos = time.time()
    print('starting integration')
    outarray = odeint(mc.phytomftm_extendedoutput_forcing, initcond, timedays_model, args=(all_params, p, z, fx))#, rtol=1e-12, atol=1e-12)
    tos1 = time.time()
    print('finished after %4.3f sec' % (tos1 - tos))

    return outarray


def callmodelrun(pfn,zn, forcing):
    # number of phytoplankton func types
    standardparams.add('pfun_num', value=pfn, vary=False)
    # number of zooplankton groups
    standardparams.add('zoo_num', value=zn, vary=False)

    if pfn == 4 and zn == 2:
         print('4P2Z - prelim model')
         all_params = (standardparams + ptype1 + ptype2 + ptype3 + ptype4 + ztype1 + ztype2)
    elif pfn == 5 and zn == 2:
         print('4P2Z - prelim model')
         all_params = (standardparams + ptype1 + ptype2 + ptype3 + ptype4 + ptype5 + ztype1 + ztype2)
    elif pfn == 2 and zn == 2:
         print('2P2Z')
         all_params = (standardparams + ptype1 + ptype2 + ztype1 + ztype2)
    elif pfn == 2 and zn == 1:
         print('2P1Z')
         all_params = (standardparams + ptype1 + ptype2 + ztype1)
    elif pfn == 1 and zn == 2:
         print('2P1Z')
         all_params = (standardparams + ptype1 + ztype1 + ztype2)
    elif pfn == 1 and zn == 1:
         print('1P1Z')
         all_params = (standardparams + ptype1 + ztype1)
    else:
        print('just standard params')
        all_params = (standardparams)

    #all_params = (standardparams)
    parameters = all_params
    initialcond = setupinitcond(pfn,zn)
    print(initialcond)
    out = runmodel(parameters,initialcond, forcing)

    return out



timedays_model = np.arange(0., 5 * 365., 1.0)

out5P2Z = callmodelrun(5,2, 'variableMLD')

out5P2Z_2 = callmodelrun(5,2,'varMLDconstNuts')

out5P2Z_3 = callmodelrun(5,2,'constantMLD')