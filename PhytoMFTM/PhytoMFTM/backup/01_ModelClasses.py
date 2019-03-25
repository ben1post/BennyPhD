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

        self.grazepref = [1 for i in range(self.zn)]
        #self.grazepref = [1 /(self.pfn*self.zn) for i in range(self.zn)]
        #self.grazepref = [1 /(self.zn) for i in range(self.pfn)]
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

    def gains(self, nuptake, siuptake, lighthrv, tempdepgro):
        if self.U_Si == 0:
            # non-diatoms
            Gains = (self.muP * nuptake * lighthrv * tempdepgro)
        else:
            # diatoms
            Gains = (self.muP * min(nuptake, siuptake) * lighthrv * tempdepgro)
        return Gains

    def grazedphyto(self, Itots, P, RperZ):
        Ri = [self.grazepref[i] * P for i in range(self.zn)]
        GrzLoss = sum([Itots[i] * Ri[i] / RperZ[i] for i in range(self.zn)])
        return GrzLoss

    def sinking(self, intMLD):
        Sink = self.v / intMLD  # Phytoplankton sinking as a function of MLD and sinking rate
        return Sink

    def grazed(self, grzing):
        Grazed = grzing
        return Grazed

    def mixing(self, diffmix):
        PMix = diffmix
        return PMix

    def silicatedrawdown(self, Gainz):
        if self.U_Si == 0:
            # non-diatoms
            return 0
        else:
            # diatoms
            Si_Drawdown = self.ratioSi * Gainz
            return Si_Drawdown

    def mortality(self,):
        Mortal = self.moP
        return Mortal


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

        self.feedpref = [1 for i in range(self.pfn)]
        #self.feedpref = [1/(self.pfn) for i in range(self.pfn)]
        #self.feedpref = [1/(self.pfn*self.zn) for i in range(self.pfn)]
        print('feedpref',self.feedpref)

    def zoomortality(self, Z):
        # Quadratic loss term for closure
        return self.moZ * Z ** 2

    def ressourcedensity(self, P):
        R = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])  # add up total phy available for this ztype
        return R

    def zoograzing(self, RperZs, P):
        Ri = [self.feedpref[i] * P[i] for i in range(self.pfn)]  # calculate available from ind phy types
        Itot = sum([Ri[i]/RperZs * RperZs / (self.Kp + RperZs) * self.muZ for i in range(self.pfn)])
        return Itot

    def zoogrowth(self, Itots):
        Growth = self.deltaZ * Itots
        return Growth

    def unassimilatedfeeding(self, Itots):
        UnAsFeeding = (1. - self.deltaZ) * Itots
        return UnAsFeeding