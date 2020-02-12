#!/usr/bin/env python
# -*- coding: utf-8 -*-

import time
import numpy as np
from lmfit import Parameters
from phydra.model import cariaco
from phydra.core import ModelSetup

from scipy.integrate import odeint

# TODO:
#  - import necessary functions & objects

######### PARAMETER SETUP #############

parameters = Parameters()

# NUTRIENT(s)
parameters.add('nuts_num', value=1)
parameters.add('nuts1_nuttype', value=0)  # Nitrate

# PHYTOPLANKTON(s)
parameters.add('phyto_num', value=5)

parameters.add('OptI', value=40, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)
#parameters.add('alpha', value=0.15, vary=False)  # initial slope of the P-I curve 0.034
#parameters.add('VpMax', value=2.5, vary=False)    # maximum photosynthetic rate

parameters.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)


#parameters.add('ratioSi', value=0, vary=False)  # Silicate ratio ## FALLBACK PARAM FOR OTHER PFTs
parameters.add('U_N', value=1.5, vary=False)    # Nitrate Half Saturation Constant
#parameters.add('U_P', value=0, vary=False)    # Phosphate Half Saturation Constant
#parameters.add('U_Si', value=0, vary=False)   # Silicate Half Saturation Constant
parameters.add('muP', value=1.2, vary=False)    # Phytoplankton maximum growth rate (d^-1)
parameters.add('v', value=1, vary=False)    # Phytoplankton sinking rate (m d^-1)

# PFTs:
# DIATOMS
parameters.add('phyto1_muP', value=1.25, vary=False)
parameters.add('phyto1_U_N', value=1.78, vary=False)
parameters.add('phyto1_OptI', value=120/3, vary=False)
parameters.add('phyto1_v', value=1.5, vary=False)

# HAPTO
parameters.add('phyto2_muP', value=0.9, vary=False)
parameters.add('phyto2_U_N', value=0.15, vary=False)
parameters.add('phyto2_OptI', value=125/3, vary=False)

# CYANO
parameters.add('phyto3_muP', value=1.4, vary=False)
parameters.add('phyto3_U_N', value=1.15, vary=False)
parameters.add('phyto3_OptI', value=300/3, vary=False)

# DINO
parameters.add('phyto4_muP', value=0.5, vary=False)
parameters.add('phyto4_U_N', value=3.15, vary=False)
parameters.add('phyto4_OptI', value=200/3, vary=False)

# OTHERS
parameters.add('phyto5_muP', value=0.7, vary=False)
parameters.add('phyto5_U_N', value=1.15, vary=False)
parameters.add('phyto5_OptI', value=100/3, vary=False)

# ZOOPLANKTON(s)
parameters.add('zoo_num', value=1)
parameters.add('moZ', value=0.01, vary=False)        # Zooplankton mortality (d^-1)
parameters.add('deltaZ', value=0.75, vary=False)    # Zooplankton Grazing assimilation coefficient (-)

parameters.add('Kp', value=1.5, vary=False)     # Zooplankton Grazing saturation constant (-)
parameters.add('pred', value=0.1, vary=False)  # quadratic higher order predation rate on zooplankton
parameters.add('muZ', value=1.0, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# ZOO Feed Prefs
parameters.add('zoo1_P1', value=.34, vary=False)
parameters.add('zoo1_P2', value=.33, vary=False)
parameters.add('zoo1_P3', value=.33, vary=False)
parameters.add('zoo1_P4', value=.33, vary=False)
parameters.add('zoo1_P5', value=.33, vary=False)
parameters.add('zoo1_D1', value=.33, vary=False)

# CONVERT FEEDPREFS TO GRAZEPREF FOR CALCULATION OF GRAZING
#parameters.add('zoo1_Zint_grazed1', value=parameters['zoo1_Zint_feed1'].value, vary=False)

parameters.add('phyto1_Z1', value=parameters['zoo1_P1'].value, vary=False)
parameters.add('phyto2_Z1', value=parameters['zoo1_P2'].value, vary=False)
parameters.add('phyto3_Z1', value=parameters['zoo1_P3'].value, vary=False)
parameters.add('phyto4_Z1', value=parameters['zoo1_P4'].value, vary=False)
parameters.add('phyto5_Z1', value=parameters['zoo1_P5'].value, vary=False)
parameters.add('det1_Z1', value=parameters['zoo1_D1'].value, vary=False)

# DETRITUS (non-biotic pools)
parameters.add('det_num', value=1)
parameters.add('deltaD_N', value=0.01, vary=False)   # Nitrate Remineralization rate (d^-1)

# PHYSICS
#parameters.add('kappa', value=0.1, vary=False)  # vary=False) # min=0.09, max=0.11) # Diffusive mixing across thermocline (m*d^-1)
parameters.add('kw', value=0.04, vary=False)     # Light attenuation constant of water (m^-1)
parameters.add('kc', value=0.03, vary=False)      # Light attenuation via phytoplankton pigment (m^-1)

#NEW EMPOWER:
#parameters.add('moP_quad', value=0.025, vary=False)    # Phytoplankton mortality (d^-1)

parameters.add('wmix', value=0.1, vary=False)
parameters.add('beta_feed', value=0.69, vary=False)
parameters.add('kN_feed', value=0.75, vary=False)
parameters.add('vD', value=1., vary=False)


######### MODEL EVALUATION CODE #############

ms = ModelSetup(parameters, physics='Box', forcing='aggTS', time='regime1')

n, p, z, d = ms.classes

physx = ms.physics

N0 = 5.
P0 = 0.1
Z0 = 0.1
D0 = 0.1

initnut = [N0 for i in range(n.num)]
initphy = [P0 for i in range(p.num)]
initzoo = [Z0 for i in range(z.num)]
initdet = [D0 for i in range(d.num)]
initout = [0. for i in range(23)]
initcond = np.concatenate([initnut, initphy, initzoo, initdet,initout], axis=None)

print(type(initcond))
timedays = np.arange(0., 5 * 365., 1.0)

# INTEGRATE:
tos = time.time()
print('starting integration')
outarray = odeint(cariaco, initcond, timedays, args=(ms, None))  # for some reason need to pass 2 args
tos1 = time.time()
print('finished after %4.3f sec' % (tos1 - tos))

#print(outarray)
