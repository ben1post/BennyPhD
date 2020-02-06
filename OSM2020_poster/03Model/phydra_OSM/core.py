#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np


# TODO:
#  - put sliceparams and similar functions into core.py, import here.
#  - dynamically create list of state variables
#  - add to forcing.py, import necessary functions here.
#  -



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
        """ dynamically create list here!"""
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
    """ This can be defined in core, but draws from forcing.py"""
    def __init__(self,params,fxtype):
        self.parameters = params
        self.type = fxtype
        if self.type == 'slab':
            self.forcing = Forcing('WOA2018', 47, -20, 1.5)
        elif self.type == 'EMPOWER':
            self.forcing = Forcing('EMPOWER')

    def K(self, MLD, mix='h+'):
        if mix == 'h+':
            MLDt = MLD[0]
            MLDt_deriv = MLD[1]
            kappa = self.parameters['kappa'].value
            return (kappa + max(MLDt_deriv, 0)) / MLDt
        if mix == 'Z':
            MLDt = MLD[0]
            MLDt_deriv = MLD[1]
            return MLDt_deriv / MLDt

    def omegaMix(self, MLD, type='std'):
        if type == 'std':
            return (self.parameters['wmix'].value + max(MLD[1],0)) / MLD[0]
        elif type == 'D':
            return (self.parameters['wmix'].value + max(MLD[1],0) + self.parameters['vD'].value) / MLD[0]

    def MLD(self,t):
        return np.array([self.forcing.MLD.return_interpvalattime(t), self.forcing.MLD.return_derivattime(t)])

    def N0(self, t):
        return self.forcing.NOX.return_interpvalattime(t)

    def PAR(self, t):
        return self.forcing.PAR.return_interpvalattime(t)

    def Tmld(self, t):
        return self.forcing.SST.return_interpvalattime(t)


class ModelSetup:
    """this needs to be dynamically collected, right?"""
    def __init__(self, params):
        self.nutrients = StateVariables(params, 'nuts')
        self.phytoplankton = StateVariables(params, 'phyto')
        self.zooplankton = StateVariables(params, 'zoo')
        self.detritus = StateVariables(params, 'det')
        self._classes = [self.nutrients, self.phytoplankton, self.zooplankton, self.detritus]

        self.physics = Physics(params, 'EMPOWER')  # 'slab' as fxtype instead

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
        rest = x[n+p+z+d:]
        return nuts, phyto, zoo, det, rest