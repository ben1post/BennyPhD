#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas
import numpy as np
import scipy.interpolate as intrp
from PhytoMFTM.AuxFuncs import sliceparams, sliceoffparams, checkreplaceparam
from PhytoMFTM.ModelClasses import Forcing
from lmfit import minimize, Parameters, Parameter, report_fit

from scipy.io import netcdf
import os

class Nutrient:
    def __init__(self, allpars, slicedpars):
        self.nuttype = checkreplaceparam(allpars, slicedpars, 'nuttype')
        self.type = self.returnnuttype()

    def mixing(self, N0, N, K):
        return K * (N0 - N)

    def returnnuttype(self):
        if self.nuttype.value == 0:
            return 'Nitrate'
        if self.nuttype.value == 1:
            return 'Phosphate'
        if self.nuttype.value == 2:
            return 'Silicate'
        if self.nuttype.value == 3:
            return 'Iron'

class Phytoplankton:
    def __init__(self, allpars, slicedpars):
        pass

    def uptake(self, N):
        pass

    def lightharvesting(self, PAR):
        pass

    def tempdepgrowth(self, Tmld):
        pass

    def mortality(self, P):
        pass

    def zoograzing(self, Gj, P, Z):
        pass

    def sinking(self, MLD, P):
        pass

    pass

class Zooplankton:
    def __init__(self, allpars, slicedpars):
        pass

    def zoofeeding(self, P, Z, func='anderson'):
        pass

    def fullgrazing(self, Gj, P, Z):
        pass
    pass

    def assimgrazing(self, ZooFeeding):
        pass

    def mortality(self, Z):
        pass


class Detritus:
    def __init__(self, allpars, slicedpars):
        pass

    def remineralisation(self, D):
        pass

    pass


class StateVariables:
    """"""
    def __init__(self,params, SVtype):
        self.type = SVtype
        self.num = params[SVtype + '_num'].value
        self.allpars = params
        self.svs = self.createlistofsvs()

    def __getattr__(self, key):
        """ This function is necessary for the StateVariables class to
        pass functions to the contained state variables when called within the ODE"""
        def fn(*args,**kwargs):
            return np.concatenate([getattr(x, key)(*args,**kwargs) for x in self.svs])
        #print(fn)
        return fn

    def sv(self, *args):
        if self.type == 'nuts':
            return Nutrient(*args)
        elif self.type == 'phyto':
            return Phytoplankton(*args)
        elif self.type == 'zoo':
            return Zooplankton(*args)
        elif self.type == 'det':
            return Detritus(*args)

    def createlistofsvs(self):
        return [self.sv(self.allpars, sliceparams(self.allpars, self.type + str(i + 1))) for i in range(self.num)]


class Physics:
    """"""
    def __init__(self,type):
        self.type = type
        self.forcing = Forcing('WOA2018',47,-20,1.5)

    def K(self,t, mix='h+'):
        if mix=='h+':
            return self.forcing.MLD.return_interpvalattime(t)
        if mix=='Z':
            return self.forcing


    def MLD(self,t):

        return [self.forcing.MLD.return_interpvalattime(t), self.forcing.MLD.return_derivattime(t)]

    def N0(self,t):
        pass

    def PAR(self,t):
        pass

    def Tmld(self,t):
        pass
    pass


class ModelSetup:
    """"""
    def __init__(self, params):
        self.nutrients = StateVariables(params, 'nuts')
        self.phytoplankton = StateVariables(params, 'phyto')
        self.zooplankton = StateVariables(params, 'zoo')
        self.detritus = StateVariables(params, 'det')
        self._classes = [self.nutrients, self.phytoplankton, self.zooplankton, self.detritus]

        self.physics = Physics(params, 'slab')

    @property
    def classes(self):
        return self._classes

    def timestep_init(self,x):
        n = self.nutrients.num
        p = self.phytoplankton.num
        z = self.zooplankton.num
        d = self.detritus.num
        nuts = x[0:n]
        phyto = x[n:n+p]
        zoo = x[n+p:n+p+z]
        det = x[n+p+z:n+p+z+d]
        return [nuts, phyto, zoo, det]


# MODEL ODE

def ode(x,t,modelsetup):
    """System of ODEs"""
    # instead of passing params, to initialize zn and pn, this should be done wrapped in a class
    # modelsetup can be creating by adding all relevant info to params object, and calling ms_init(params)

    N,P,Z,D = modelsetup.timestep_init(x)

    # N = [Ni,P,Si]
    # P = [P1,P2,P3]
    # Z = [Z1,Z2,Z3]
    # D = [D]

    # 'n' contains all functions related to the different nutrients
    # a call to a function goes through a np.array of all nuts with that specific function,
    # always using the respective param set, and returning the results as that same np.array.

    MLD = physx.MLD(t) # MLD = [int_MLD, deriv_MLD]
    N0 = physx.N0(t) # N0 = [Ni0,P0,Si0]
    PAR = physx.PAR(t)
    Tmld = physx.Tmld(t)

    # the functions have to return an np.array!! this way piecewise adding is possible! no list
    K = physx.K(MLD) # i.e. there is constant mixing & increased mix when MLD shallowing
    K_Z = physx.K(MLD, mix='Z')

    # Grazing
    Gj = z.zoofeeding(P, Z, func='anderson')  # feeding probability for all food
    ZooFeeding = z.fullgrazing(Gj, P, Z)

    # Phytoplankton Fluxes
    PGains = p.uptake(N) * p.lightharvesting(PAR) * p.tempdepgrowth(Tmld)
    PMortality = p.mortality(P)
    PhytoLosses = p.zoograzing(Gj, P, Z) + PMortality + p.sinking(MLD, P) + P * K

    # Zooplankton Fluxes
    ZGains = z.assimgrazing(ZooFeeding)
    ZMortality = z.mortality(Z)
    ZLosses = ZMortality + z.higherorderpred(Z) + Z * K_Z

    # Detritus Fluxes
    DGains = sum(z.unassimilatedgrazing(ZooFeeding)) + sum(ZMortality) + sum(PMortality)
    DRemin = d.remineralisation(D)
    DLosses = DRemin - D * K

    # THESE ALL HAVE TO BE np.arrays, otherwise this won't work!

    # ALSO using SUM() will cause problems with the flexibility!!!

    Px = PGains - PhytoLosses

    Nx = DRemin + n.mixing(N0, N, K) - sum(PGains)  # Nutrient draw down

    Zx = ZGains - ZLosses  # Zooplankton losses due to mortality and mixing

    Dx = DGains - DLosses   # Detritus

    out = np.concatenate(Nx, Px, Zx, Dx)
    return out



######### PARAMETER SETUP #############

parameters = Parameters()

# NUTRIENT(s)
parameters.add('nuts_num', value=1)
parameters.add('nuts1_nuttype', value='0') # Nitrate

# PHYTOPLANKTON(s)
parameters.add('v', value=0.01, vary=False)      # Sinking of Phytoplankton from Mixed Layer
parameters.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)
parameters.add('alpha', value=0.15, vary=False)  # initial slope of the P-I curve
parameters.add('VpMax', value=1., vary=False)    # maximum photosynthetic rate

parameters.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)

parameters.add('ratioSi', value=0, vary=False)  # Silicate ratio ## FALLBACK PARAM FOR OTHER PFTs
parameters.add('U_N', value=0, vary=False)    # Nitrate Half Saturation Constant
parameters.add('U_Si', value=0, vary=False)   # Silicate Half Saturation Constant
parameters.add('muP', value=0, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# ZOOPLANKTON(s)
parameters.add('moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
parameters.add('deltaZ', value=0.75, vary=False)    # Zooplankton Grazing assimilation coefficient (-)

parameters.add('Kp', value=0, vary=False)     # Zooplankton Grazing saturation constant (-)
parameters.add('pred', value=0, vary=False)  # quadratic higher order predation rate on zooplankton
parameters.add('muZ', value=0, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# DETRITUS (non-biotic pools)
parameters.add('deltaD_N', value=0.1, vary=False)   # Nitrate Remineralization rate (d^-1)

# PHYSICS
parameters.add('kappa', value=0.1, vary=False) # vary=False) # min=0.09, max=0.11) # Diffusive mixing across thermocline (m*d^-1)
parameters.add('kw', value=0.04, vary=False)     # Light attenuation constant of water (m^-1)
parameters.add('kc', value=0.03, vary=False)      # Light attenuation via phytoplankton pigment (m^-1)



######### MODEL EVALUATION CODE #############

ms = ModelSetup(parameters)

n,p,z,d = ms.classes

physx = ms.physics

print(n)

initnut = [N0, Si0, D0]
initzoo = [Z0 for i in range(zn)]
initphy = [P0 for i in range(pfn)]
outputl = [0 for i in range(20)]
initcond = np.concatenate([initnut, initzoo, initphy, outputl])

timedays = np.arange(0,30,1)

import time
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# INTEGRATE:
tos = time.time()
print('starting integration')
outarray = odeint(ode, initcond, timedays, args=())
tos1 = time.time()
print('finished after %4.3f sec' % (tos1 - tos))

print(outarray)

#plt.plot(outarray[0],timedays,c='r')
#plt.plot(outarray[1],timedays,c='b')