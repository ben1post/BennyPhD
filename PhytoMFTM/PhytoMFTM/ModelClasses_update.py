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
    """"""
    def __init__(self, allpars, slicedpars, num):
        self.num = num
        self.nuttype = checkreplaceparam(allpars, slicedpars, 'nuttype')
        self.type = self.returnnuttype()

    def mixing(self, N0, N, K):
        """"""
        # TODO this only works for Nitrate as of now
        return K * (N0 - N[self.num])

    def returnnuttype(self):
        if self.nuttype == 0:
            return 'Nitrate'
        if self.nuttype == 1:
            return 'Phosphate'
        if self.nuttype == 2:
            return 'Silicate'
        if self.nuttype == 3:
            return 'Iron'

class Phytoplankton:
    """"""
    def __init__(self, allpars, slicedpars, num):
        self.num = num

        self.kw = checkreplaceparam(allpars, slicedpars, 'kw')
        self.OptI = checkreplaceparam(allpars, slicedpars, 'OptI')

        self.U_N = checkreplaceparam(allpars, slicedpars, 'U_N')
        self.U_P = checkreplaceparam(allpars, slicedpars, 'U_P')
        self.U_Si = checkreplaceparam(allpars, slicedpars, 'U_Si')

        self.v = checkreplaceparam(allpars, slicedpars, 'v')

        self.muP = checkreplaceparam(allpars, slicedpars, 'muP')
        self.moP = checkreplaceparam(allpars, slicedpars, 'moP')

        self.ratioSi = checkreplaceparam(allpars, slicedpars, 'ratioSi')

        self.pfn = allpars['phyto_num'].value
        self.zn = allpars['zoo_num'].value

        self.kc = checkreplaceparam(allpars, slicedpars, 'kc')
        self.alpha = checkreplaceparam(allpars, slicedpars, 'alpha')
        self.VpMax = checkreplaceparam(allpars, slicedpars, 'VpMax')

        self.zoolist = ['Z'+str(j+1) for j in range(self.zn)]
        self.grazepref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoolist]

    def uptake(self, N):
        """Michealis-Menten uptake of Nutrients"""
        # TODO: add functionality of multiple Nutrients/formulations
        Uptake = N[0] / (N[0] + self.U_N)
        return Uptake

    def lightharvesting(self, MLD, PAR, P, type='Smith'):
        """Light - modification of phytoplankton growth rate"""
        if type == 'Steele':
            lighthrv = 1. / (self.kw * MLD) * \
                       (-np.exp(1. - PAR / self.OptI) - (
                           -np.exp((1. - (PAR * np.exp(-self.kw * MLD)) / self.OptI))))
            return lighthrv
        if type == 'Smith':
            kPAR = self.kw + self.kc * sum(P)
            x_0 = self.alpha * PAR * np.exp(- kPAR * 0)
            x_H = self.alpha * PAR * np.exp(- kPAR * MLD)
            VpH = (self.VpMax / (kPAR * MLD)) * \
                  np.log(
                      (x_0 + np.sqrt(self.VpMax ** 2 + x_0 ** 2)) /
                      (x_H + np.sqrt(self.VpMax ** 2 + x_H ** 2))
                  )
            return VpH

    def tempdepgrowth(self, Tmld):
        """"""
        tdp = np.exp(0.063 * Tmld)
        return tdp

    def mortality(self, P):
        """"""
        Mortal = self.moP * P[self.num]
        return Mortal

    def zoograzing(self, Gj, P, Z):
        """"""
        # take the general grazing term from each zooplankton, multiply by phyto fraction and sum
        Grazing = [Gj[j] * (self.grazepref[j] * P[self.num] ** 2) * Z[j] for j in range(self.zn)]
        GrazingPerZ = sum(Grazing)
        return GrazingPerZ

    def sinking(self, MLD, P):
        """"""
        Sink = self.v / MLD * P[self.num]  # Phytoplankton sinking as a function of MLD and sinking rate
        return Sink

class Zooplankton:
    """"""
    def __init__(self, allpars, slicedpars, num):
        self.num = num
        # zooplankton
        self.moZ = checkreplaceparam(allpars, slicedpars, 'moZ')
        self.muZ = checkreplaceparam(allpars, slicedpars, 'muZ')
        # grazing params
        self.Kp = checkreplaceparam(allpars, slicedpars, 'Kp')
        self.deltaZ = checkreplaceparam(allpars, slicedpars, 'deltaZ')

        #self.muIntGraze = checkreplaceparam(allpars, slicedpars, 'muIntGraze')
        #self.kIntGraze = checkreplaceparam(allpars, slicedpars, 'kIntGraze')

        self.pred = checkreplaceparam(allpars, slicedpars, 'pred')
        #self.deltaLambda = checkreplaceparam(allpars, slicedpars, 'deltaLambda')

        self.pfn = allpars['phyto_num'].value
        self.zn = allpars['zoo_num'].value

        self.beta = 2  # for Vallina KTW Grazing, beta = 1 : Holling Type II, beta = 2 : Holling Type III
        self.ksat = self.Kp

        self.phylist = ['P' + str(i + 1) for i in range(self.pfn)]
        self.zoointlistfeed = ['Zint_feed' + str(j + 1) for j in range(self.zn)]  # list of feeding
        self.zoointlistgrazed = ['Zint_grazed' + str(j + 1) for j in range(self.zn)]
        # print(self.phylist, self.zoointlistfeed, self.zoointlistgrazed)

        self.feedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.phylist]
        self.interfeedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoointlistfeed]
        self.intergrazedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoointlistgrazed]

    def zoofeeding(self, P, Z, func='anderson'):
        if func == 'anderson':
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP + FrhoZ
            GrazingProb = self.muZ / (self.ksat ** 2 + Frho)
            return GrazingProb

        elif func == 'fasham': #put the holling type num here
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP + FrhoZ  # active switching coefficient
            FpP = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])
            FpZ = sum([self.interfeedpref[j] * Z[j] for j in range(self.zn)])
            Fp = FpP + FpZ  # total food available
            GrazingProb = self.muZ * (1 / (self.ksat * Fp + Frho))
            return GrazingProb

        elif func == 'vallina':
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP + FrhoZ  # active switching coefficient
            FpP = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])
            FpZ = sum([self.interfeedpref[j] * Z[j] for j in range(self.zn)])
            Fp = FpP + FpZ  # total food available
            GrazingProb = self.muZ * (1 / Frho) * ((Fp ** self.beta) / (self.ksat ** self.beta + Fp ** self.beta))
            return GrazingProb

        else:
            print('no grazing formulation given, wrong func key')

    def fullgrazing(self, Gj, P, Z):
        # phytouptake + zooplankton per zooplankton for each phyto
        IprobP = [Gj[self.num] * (self.feedpref[i] * P[i] ** 2) for i in range(self.pfn)]  # grazeprob per each PFT
        IprobZ = [Gj[self.num] * (self.interfeedpref[j] * Z[j] ** 2) for j in range(self.zn)]
        Iprob = IprobP + IprobZ
        #perhaps here = multiply by Z before summing Iprobs!
        Itots = sum(Iprob)
        Itot = Itots * Z[self.num]
        return Itot

    def assimgrazing(self, ZooFeeding):
        AssimGrazing = self.deltaZ * ZooFeeding[self.num]
        return AssimGrazing

    def unassimilatedgrazing(self, ZooFeeding):
        UnAsGraze = (1. - self.deltaZ) * ZooFeeding[self.num]
        return UnAsGraze

    def mortality(self, Z, type='linear'):
        if type == 'linear':
            total_moZ = self.moZ * Z[self.num]
            return total_moZ
        if type == 'quadratic':
            total_moZ = self.pred * Z[self.num]
            return total_moZ


class Detritus:
    def __init__(self, allpars, slicedpars, num):
        self.num = num
        self.deltaD = checkreplaceparam(allpars, slicedpars, 'deltaD_N')

    def remineralisation(self, D):
        return D[self.num] * self.deltaD


class StateVariables:
    """"""
    def __init__(self,params, SVtype):
        self.type = SVtype
        self.num = params[SVtype + '_num'].value
        self.allpars = params
        self.svs = self.createlistofsvs()
        print(self.type,self.num,'created')

    def __getattr__(self, key):
        """ This function is necessary for the StateVariables class to
        pass functions to the contained state variables when called within the ODE"""
        def fn(*args,**kwargs):
            return np.array([getattr(x, key)(*args,**kwargs) for x in self.svs])
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
        return [self.sv(self.allpars, sliceparams(self.allpars, self.type + str(i + 1)), i) for i in range(self.num)]


class Physics:
    """"""
    def __init__(self,params,type):
        self.parameters = params
        self.type = type
        self.forcing = Forcing('WOA2018', 47, -20, 1.5)

    def K(self,t, MLD, mix='h+'):
        if mix == 'h+':
            MLDt = MLD[0]
            MLDt_deriv = MLD[1]
            kappa = self.parameters['kappa'].value
            return (kappa + max(MLDt_deriv, 0)) / MLDt
        if mix == 'Z':
            MLDt = MLD[0]
            MLDt_deriv = MLD[1]
            return MLDt_deriv / MLDt

    def MLD(self,t):
        return np.array([self.forcing.MLD.return_interpvalattime(t), self.forcing.MLD.return_derivattime(t)])

    def N0(self, t):
        return self.forcing.NOX.return_interpvalattime(t)

    def PAR(self, t):
        return self.forcing.PAR.return_interpvalattime(t)

    def Tmld(self, t):
        return self.forcing.SST.return_interpvalattime(t)


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
        return nuts, phyto, zoo, det


# MODEL ODE

def ode(x,t,modelsetup, q):
    """System of ODEs"""
    # instead of passing params, to initialize zn and pn, this should be done wrapped in a class
    # modelsetup can be creating by adding all relevant info to params object, and calling ms_init(params)

    N,P,Z,D = modelsetup.timestep_init(x)
#    print(N,P,Z,D)

    # N = [Ni,P,Si]
    # P = [P1,P2,P3]
    # Z = [Z1,Z2,Z3]
    # D = [D]

    # 'n' contains all functions related to the different nutrients
    # a call to a function goes through a np.array of all nuts with that specific function,
    # always using the respective param set, and returning the results as that same np.array.

    MLD = physx.MLD(t)  # MLD = [int_MLD, deriv_MLD]
    N0 = physx.N0(t)  # N0 = [Ni0,P0,Si0]
    PAR = physx.PAR(t)
    Tmld = physx.Tmld(t)

    # the functions have to return an np.array!! this way piecewise adding is possible! no list
    K = physx.K(t, MLD)  # i.e. there is constant mixing & increased mix when MLD shallowing
    K_Z = physx.K(t, MLD, mix='Z')

    # Grazing
    Gj = z.zoofeeding(P, Z, func='anderson')  # feeding probability for all food

    ZooFeeding = z.fullgrazing(Gj, P, Z)


    # Phytoplankton Fluxes
    PGains = p.uptake(N) * p.lightharvesting(MLD[0], PAR, P) * p.tempdepgrowth(Tmld)
    PMortality = p.mortality(P)
    PhytoLosses = p.zoograzing(Gj, P, Z) + PMortality + p.sinking(MLD[0], P) + P * K

    # Zooplankton Fluxes
    ZGains = z.assimgrazing(ZooFeeding)

    ZMortality = z.mortality(Z, type='linear')
    ZLosses = ZMortality + z.mortality(Z, type='quadratic') + Z * K_Z


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

    out = [Nx, Px, Zx, Dx]
    return np.concatenate(out)



######### PARAMETER SETUP #############

parameters = Parameters()

# NUTRIENT(s)
parameters.add('nuts_num', value=1)
parameters.add('nuts1_nuttype', value=0)  # Nitrate

# PHYTOPLANKTON(s)
parameters.add('phyto_num', value=2)
parameters.add('v', value=0.01, vary=False)      # Sinking of Phytoplankton from Mixed Layer
parameters.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)
parameters.add('alpha', value=0.15, vary=False)  # initial slope of the P-I curve
parameters.add('VpMax', value=1., vary=False)    # maximum photosynthetic rate

parameters.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)

parameters.add('ratioSi', value=0, vary=False)  # Silicate ratio ## FALLBACK PARAM FOR OTHER PFTs
parameters.add('U_N', value=0.5, vary=False)    # Nitrate Half Saturation Constant
parameters.add('U_P', value=0, vary=False)    # Phosphate Half Saturation Constant
parameters.add('U_Si', value=0, vary=False)   # Silicate Half Saturation Constant
parameters.add('muP', value=0, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# ZOOPLANKTON(s)
parameters.add('zoo_num', value=1)
parameters.add('moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
parameters.add('deltaZ', value=0.75, vary=False)    # Zooplankton Grazing assimilation coefficient (-)

parameters.add('Kp', value=1.5, vary=False)     # Zooplankton Grazing saturation constant (-)
parameters.add('pred', value=0.01, vary=False)  # quadratic higher order predation rate on zooplankton
parameters.add('muZ', value=0.5, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# DETRITUS (non-biotic pools)
parameters.add('det_num', value=1)
parameters.add('deltaD_N', value=0.05, vary=False)   # Nitrate Remineralization rate (d^-1)

# PHYSICS
parameters.add('kappa', value=0.1, vary=False)  # vary=False) # min=0.09, max=0.11) # Diffusive mixing across thermocline (m*d^-1)
parameters.add('kw', value=0.04, vary=False)     # Light attenuation constant of water (m^-1)
parameters.add('kc', value=0.03, vary=False)      # Light attenuation via phytoplankton pigment (m^-1)

# MIKRO
parameters.add('zoo1_P1', value=.7, vary=False)  # Diatoms
parameters.add('zoo1_P2', value=.3, vary=False)

# MESO

# inter zoo feeding
# MIKRO
parameters.add('zoo1_Zint_feed1', value=0, vary=False)

# CONVERT FEEDPREFS TO GRAZEPREF FOR CALCULATION OF GRAZING
parameters.add('zoo1_Zint_grazed1', value=parameters['zoo1_Zint_feed1'].value, vary=False)

parameters.add('phyto1_Z1', value=parameters['zoo1_P1'].value, vary=False)
parameters.add('phyto2_Z1', value=parameters['zoo1_P2'].value, vary=False)

######### MODEL EVALUATION CODE #############

ms = ModelSetup(parameters)

n,p,z,d = ms.classes

physx = ms.physics


N0 = 1
P0 = 0.5
Z0 = 0.1
D0 = 0.1

initnut = [N0 for i in range(n.num)]
initphy = [P0 for i in range(p.num)]
initzoo = [Z0 for i in range(z.num)]
initdet = [D0 for i in range(d.num)]
initcond = np.concatenate([initnut, initphy, initzoo, initdet], axis=None)

timedays = np.arange(0., 5 * 365., 1.0)

import time
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# INTEGRATE:
tos = time.time()
print('starting integration')
outarray = odeint(ode, initcond, timedays, args=(ms,1))
tos1 = time.time()
print('finished after %4.3f sec' % (tos1 - tos))

#print(outarray)


timedays_ly = timedays#[1:366]
# truncate outarraySiNO to last year of 5:
outarray_ly = outarray#[1460:1825]

numcols = 1
f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, numcols, sharex='col')#, sharey='row')

ax1.set_title('N')
ax1.plot(timedays_ly,outarray_ly[:,0],label='N')
ax1.set_ylim(bottom=0)

ax2.set_title('P')
ax2.plot(timedays_ly,outarray_ly[:,1],label='P1')
ax2.plot(timedays_ly,outarray_ly[:,2],label='P2')
ax2.legend()
ax2.set_ylim(bottom=0)

ax3.set_title('Z')
ax3.plot(timedays_ly,outarray_ly[:,3],label='Z')
ax3.set_ylim(bottom=0)

ax4.set_title('D')
ax4.plot(timedays_ly,outarray_ly[:,4],label='D')
ax4.set_ylim(bottom=0)

plt.tight_layout()
