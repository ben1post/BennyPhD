#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from phydra.aux import sliceparams, sliceoffparams
from phydra.classes import Nutrient, Phytoplankton, Zooplankton, Detritus
from phydra.forcing import Forcing

# TODO:
#  - dynamically create list of state variables
#  -

class StateVariables:
    """Composite class that handles all instances of classes in phydra.classes
    contains list of instances of respective class, and handles calling function on each instance,
    when called within system of ODEs
    """
    # TODO: AS OF NOW PFT params need to be situatied in dict AFTER stdparams, this needs to be fixed
    def __init__(self, params, SVtype):
        self.type = SVtype
        self.num = params[SVtype + '_num'].value

        self.allpars = params

        self.svs = self.createlistofsvs()


    def __getattr__(self, key):
        """This function is necessary for the StateVariables class to
        pass functions to the contained state variable(s) when called within the ODE"""
        def fn(*args,**kwargs):
            return np.array([getattr(x, key)(*args, **kwargs) for x in self.svs])
        return fn

    def sv(self, *args):
        """Called by createlistofsvs function below, handles creating correct instance based on parameters passed"""
        if self.type == 'nuts':
            return Nutrient(*args)
        elif self.type == 'phyto':
            return Phytoplankton(*args)
        elif self.type == 'zoo':
            return Zooplankton(*args)
        elif self.type == 'det':
            return Detritus(*args)

    def createlistofsvs(self):
        """Function to create list of class instances, e.g. 2 Phytoplankton class instances"""
        return [self.sv(self.allpars, sliceparams(self.allpars, self.type + str(i + 1)), i) for i in range(self.num)]


class Physics:
    """Composite class handles model forcing, contains forcing and verification data,
    actual forcing & verification is defined in phydra.forcing module"""
    def __init__(self,params, phsxtype, fxtype, time, pad, extendstart):
        self.parameters = params
        self.type = phsxtype
        if self.type == 'EMPOWER':
            self.forcing = Forcing('EMPOWER')
        elif self.type == 'Box':
            self.forcing = Forcing(fxtype, time, pad, extendstart)
            self.BoxDepth = 100

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

    def wMix(self, X258, type='std'):
        if type == 'std':
            return (self.parameters['wmix'].value + max(-X258[1], 0)) / self.BoxDepth

    def MLD(self,t):
        return np.array([self.forcing.MLD.return_interpvalattime(t), self.forcing.MLD.return_derivattime(t)])

    def X258(self, t):
        return np.array([self.forcing.X258.return_interpvalattime(t), self.forcing.X258.return_derivattime(t)])

    def N0(self, t):
        # TODO: figure out a better way to deal with nutrients (add more!?)
        return [self.forcing.NOX.return_interpvalattime(t)]

    def Si0(self, t):
        # TODO: figure out a better way to deal with nutrients (add more!?)
        return [self.forcing.SiOH.return_interpvalattime(t)]

    def PAR(self, t):
        return self.forcing.PAR.return_interpvalattime(t)

    def Tmld(self, t):
        return self.forcing.SST.return_interpvalattime(t)


class ModelSetup:
    """Composite class containing all relevant classes related to model setup,
    can be passed to ode function and to plotting scripts to access state variables, forcing and verification data"""

    def __init__(self, params, physics='Box', forcing='aggTS', time='regime1', pad=False, extendstart=False):
        self.nutrients = StateVariables(params, 'nuts')
        self.phytoplankton = StateVariables(params, 'phyto')
        self.zooplankton = StateVariables(params, 'zoo')
        self.detritus = StateVariables(params, 'det')
        self._classes = [self.nutrients, self.phytoplankton, self.zooplankton, self.detritus]

        self.physics = Physics(params, physics, forcing, time, pad, extendstart)  # 'slab' as fxtype instead

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
        _ = x[n+p+z+d:]
        return nuts, phyto, zoo, det, _