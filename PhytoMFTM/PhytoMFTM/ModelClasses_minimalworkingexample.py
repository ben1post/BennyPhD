#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas
import numpy as np
import scipy.interpolate as intrp
from PhytoMFTM.AuxFuncs import sliceparams, sliceoffparams, checkreplaceparam

from lmfit import minimize, Parameters, Parameter, report_fit

from scipy.io import netcdf
import os

class Nutrient:
    def __init__(self, allpars, slicedpars):
        self.kappa = checkreplaceparam(allpars, slicedpars, 'kappa')

    def mixing(self, N):
        return self.kappa * N

"""
class SV:
    def __init__(self, params, SVStype, SVparams):
        #Class that contains all general properties of state variables, to be inherited by each instance
        self.type = SVStype
        self.sliced_pars = SVparams
        self.all_pars = params
        if SVStype == 'nuts':
            self.funcs = Nutrient(self.all_pars,self.sliced_pars)
        #if SVStype == 'phyto':
        #    self.funcs = Phyto(self.all_pars,self.pars)
"""

class StateVariables:
    """"""
    def __init__(self,params, SVtype):
        self.type = SVtype
        self.num = params[SVtype + '_num'].value
        self.allpars = params
        self.svs = self.createlistofsvs()

    def __getattr__(self, key):
        def fn(*args,**kwargs):
            return np.concatenate([getattr(x, key)(*args,**kwargs) for x in self.svs])
        #print(fn)
        return fn

    def createlistofsvs(self):
        if self.type == 'nuts':
            print([self.type + str(i + 1) for i in range(self.num)])
            return [Nutrient(self.allpars, sliceparams(self.allpars, self.type + str(i + 1))) for i in range(self.num)]

    def init(self):
        return self.svs

class ModelSetup:
    def __init__(self, params):
        self.nutrients = StateVariables(params, 'nuts')
        #self.phytoplankton = Plankton('in here unpack the number of respective instances, return list of instances')
        self._classes = self.nutrients  # concatenate instances here ,self.phytoplankton,self.zooplankton,self.detritus

    @property
    def classes(self):
        return self._classes




def ode(x,t):
    N = x

    Ngain = n.mixing(N)

    #print(Ngain)
    Nx = sum(Ngain)

    return Nx


parameters = Parameters()

parameters.add('nuts_num', value=2)

parameters.add('kappa', value=0.)

parameters.add('nuts1_kappa', value=0.1)
parameters.add('nuts10_kappa', value=0.1)

ms = ModelSetup(parameters)

n = ms.classes
print(n)

initcond = 1

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