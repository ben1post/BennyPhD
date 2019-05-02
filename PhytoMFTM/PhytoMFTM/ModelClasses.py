#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from PhytoMFTM.AuxFuncs import sliceparams,sliceoffparams,checkreplaceparam

class Plankton:
    """
    initializes the Plankton (both Zoo- & Phyto-) community according to number of types prescribed by params
    """

    def __init__(self, modelparams, zptype):
        if zptype == 'Zooplankton':
            self.zoonum = modelparams['zoo_num'].value
            self.stdparams = sliceoffparams(modelparams, 'zt')
            self.pars = [ZooType(self.stdparams, sliceparams(modelparams, 'zt' + str(i + 1)))
                        for i in range(self.zoonum)]
        elif zptype == 'Phytoplankton':
            self.phytonum = modelparams['pfun_num'].value
            self.stdparams = sliceoffparams(modelparams, 'pt')
            self.pars = [PhytoType(self.stdparams, sliceparams(modelparams, 'pt' + str(i + 1))) for i in
                         range(self.phytonum)]

        else:
            raise('wrong functional type passed to Plankton class')

    def init(self, ):
        return self.pars

class PhytoType:
    """sets up individual functional types and defines functions
    such as nutrient uptake and growth"""

    def __init__(self, stdpars, slicedpars):
        self.kw = checkreplaceparam(stdpars, slicedpars, 'kw')
        self.OptI = checkreplaceparam(stdpars, slicedpars, 'OptI')

        self.U_N = checkreplaceparam(stdpars, slicedpars, 'U_N')
        self.U_Si = checkreplaceparam(stdpars, slicedpars, 'U_Si')

        self.v = checkreplaceparam(stdpars, slicedpars, 'v')

        self.muP = checkreplaceparam(stdpars, slicedpars, 'muP')
        self.moP = checkreplaceparam(stdpars, slicedpars, 'moP')

        self.ratioSi = checkreplaceparam(stdpars, slicedpars, 'ratioSi')

        self.pfn = stdpars['pfun_num']
        self.zn = stdpars['zoo_num']

        self.grazepref = [checkreplaceparam(stdpars, slicedpars, string) for string in ['Z1','Z2']]
        #self.grazepref = [1 /(self.pfn*self.zn) for i in range(self.zn)]
        #self.grazepref = [1 /(self.zn) for i in range(self.zn)]
        print('grazepref', self.grazepref)


    def n_uptake(self, Nitrate):
        N_Uptake = Nitrate / (Nitrate + self.U_N)  # Michaelis Menten - uptake of Nitrate
        return N_Uptake

    def si_uptake(self, Silicate):
        if self.U_Si == 0:
            # non-diatoms
            return 0
        else:
            # diatoms
            Si_Uptake = Silicate / (Silicate + self.U_Si)  # Michaelis Menten - uptake of Nitrate
            return Si_Uptake

    def lightharvesting(self, intMLD, intPAR):
        lighthrv =  1. / (self.kw * intMLD) * \
               (-np.exp(1. - intPAR / self.OptI) - (-np.exp((1. - (intPAR * np.exp(-self.kw * intMLD)) / self.OptI))))
        return lighthrv

    def tempdepgrowth(self, int_SST):
        tdp = np.exp(0.063 * int_SST)
        return tdp

    def gains(self, nuptake, siuptake, lighthrv, tempdepgro, P):
        if self.U_Si == 0:
            # non-diatoms
            Gains = (self.muP * nuptake * lighthrv * tempdepgro) * P
        else:
            # diatoms
            Gains = (self.muP * min(nuptake, siuptake) * lighthrv * tempdepgro) * P
        return Gains

    def sinking(self, intMLD, P):
        Sink = self.v / intMLD * P  # Phytoplankton sinking as a function of MLD and sinking rate
        return Sink

    def grazed(self, grzing, P):
        Grazed = grzing #* P
        return Grazed

    def mixing(self, diffmix, P):
        PMix = P * diffmix
        return PMix

    def silicatedrawdown(self, Gainz):
        if self.U_Si == 0:
            # non-diatoms
            return 0
        else:
            # diatoms
            Si_Drawdown = self.ratioSi * Gainz
            return Si_Drawdown

    def mortality(self, P):
        Mortal = self.moP * P
        return Mortal

    def zoograzing(self, Itot, Rj, Pi, Z):
        Grazing = sum([Itot[j] * (self.grazepref[j] * Pi) / Rj[j] * Z[j] for j in range(self.zn)])
        return Grazing


class ZooType:
    def __init__(self, stdpars, slicedpars):
        # zooplankton
        self.moZ = checkreplaceparam(stdpars, slicedpars, 'moZ')
        self.muZ = checkreplaceparam(stdpars, slicedpars, 'muZ')
        # grazing params
        self.gr_p = checkreplaceparam(stdpars, slicedpars, 'gr_p')
        self.Kp = checkreplaceparam(stdpars, slicedpars, 'Kp')
        self.deltaZ = checkreplaceparam(stdpars, slicedpars, 'deltaZ')

        self.pfn = stdpars['pfun_num']
        self.zn = stdpars['zoo_num']
        #self.muZ = self.muZ / self.zn


        self.feedpref = [checkreplaceparam(stdpars, slicedpars, string) for string in ['P1','P2','P3','P4']]
        #self.feedpref = [1 for i in range(self.pfn)]

        print('feedpref',self.feedpref, self.muZ)

    def zoomortality(self, Z):
        # Quadratic loss term for closure --> now linear, quadratic loss via higher order pred (below)
        # The problem was here! now calc as aggregate mortality and split up between types
        total_moZ = self.moZ/10 * Z #** 2  # sum(Z)
        return total_moZ#/self.zn

    def ressourcedensity(self, P):
        R = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])  # add up total phy available for this ztype
        return R

    def itot(self, Rtoti):
        Itoti = self.muZ * Rtoti/(Rtoti+self.Kp)  # add up total phy available for this ztype
        return Itoti

    def assimgrazing(self, Itot, Z):
        AssimGrazing = self.deltaZ * Itot * Z
        return AssimGrazing

    def unassimilatedgrazing(self, Itots, Z):
        UnAsGraze = (1. - self.deltaZ) * Itots * Z
        return UnAsGraze

    def interzoograze(self,i, Z):
        totalgraze= Z[0] * 0.1 * Z[1]
        if i == 0:
            return -totalgraze
        if i == 1:
            return totalgraze
        else:
            return 0

    def higherorderpred(self,Zi):
        total_moZ = self.moZ * Zi ** 2  # sum(Z)
        return total_moZ / self.zn
