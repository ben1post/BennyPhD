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

        self.zn = stdpars['zoo_num']

        self.grazepref = [1 for i in range(self.zn)]



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
        return 1. / (self.kw * intMLD) * \
               (-np.exp(1. - intPAR / self.OptI) - \
                (-np.exp((1. - (intPAR * np.exp(-self.kw * intMLD)) / self.OptI))))

    def gains(self, nuptake, siuptake, lighthrv, tempdepgro):
        if self.U_Si == 0:
            # non-diatoms
            Gains = self.muP * nuptake * lighthrv * tempdepgro
        else:
            # diatoms
            Gains = self.muP * min(nuptake, siuptake) * lighthrv * tempdepgro
        return Gains



    def grazedphyto(self, Itot, P, R):   #NEEEWW GRAZED PHYTO FORMULATION USING GRAZEDPREF LIST

        Ri = [self.grazepref[i]*P for i in range(self.zn)] # this calculates the resource density for this phyto type as grazed by all zoo types
        # Ri[1] is the relative contribution of this phyto type to Z type 1 's diet
        GrzLoss = sum([Itot[i] * (Ri[i] / R[i]) for i in range(self.zn)]) # the Itot relates to each Zootype, similarly the R should relate to each Zoo type

        return GrzLoss



    def losses(self, intMLD, grzing, diffmix):
        Sinking = self.v / intMLD  # Phytoplankton sinking as a function of MLD and sinking rate
        OtherPMortalities = self.moP  # Linear Phytoplankton mortality
        Losses = grzing + Sinking + OtherPMortalities + diffmix
        return Losses

    def silicatedrawdown(self, P, Gainz):
        Si_Drawdown = P * self.ratioSi * Gainz
        return Si_Drawdown

    def mortality(self, P):
        Mortal = self.moP * P
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

        self.feedpref = [1 for i in range(self.pfn)]

    def zoomortality(self, Z):
        # Quadratic loss term for closure
        return self.moZ * Z ** 2

    def zoograzing(self, P):  #NEWWWW formulation incl. feedpref list
        R = sum([self.feedpref[i] * P[i] for i in range(self.pfn)]) # this should add together the total ressource density per zoo type

        Itot = (R / (self.Kp + R)) * self.muZ  # this should return the
        return Itot, R

    def zoogrowth(self, Itots):
        Growth = self.deltaZ * Itots
        return Growth

    def unassimilatedfeeding(self, Itots):
        UnAsFeeding = (1. - self.deltaZ) * Itots
        return UnAsFeeding

# add grazing preference parameter
