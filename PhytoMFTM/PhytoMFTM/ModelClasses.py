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
        print(self.pfn, self.zn)

        self.kc = checkreplaceparam(stdpars, slicedpars, 'kc')
        self.alpha = checkreplaceparam(stdpars, slicedpars, 'alpha')
        self.VpMax = checkreplaceparam(stdpars, slicedpars, 'VpMax')

        self.zoolist = ['Z'+str(j+1) for j in range(self.zn)]
        print(self.zoolist)
        self.grazepref = [checkreplaceparam(stdpars, slicedpars, string) for string in self.zoolist]

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
        # light harvesting according to Steele PI curve, and a constant extinction coefficient
        self.kw = 0.1
        lighthrv =  1. / (self.kw * intMLD) * \
               (-np.exp(1. - intPAR / self.OptI) - (-np.exp((1. - (intPAR * np.exp(-self.kw * intMLD)) / self.OptI))))
        return lighthrv

    def smithpi(self, intMLD, intPAR, P):
        # lightharvesting according to the smith equation (from EMPOWER publication)
        # still need to add kPAR calculation based on kw and kchla
        kPAR = self.kw + self.kc * sum(P)
        x_0 = self.alpha * intPAR * np.exp(- kPAR * 0)
        x_H = self.alpha * intPAR * np.exp(- kPAR * intMLD)
        VpH = (self.VpMax / (kPAR * intMLD)) * \
              np.log(
                  (x_0 + np.sqrt(self.VpMax ** 2 + x_0 ** 2)) /
                  (x_H + np.sqrt(self.VpMax ** 2 + x_H ** 2))
              )
        return VpH

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

    def zoograzing(self, Gj, Pi, Z, j):
        # take the general grazing term from each zooplankton, multiply by phyto fraction and sum
        Grazing = [Gj[j] * (self.grazepref[j] * Pi ** 2) * Z[j] for j in range(self.zn)]
        GrazingPerZ = sum(Grazing)

        print('p1', Pi, 'p2', j, 'Z', Z, 'Grazing', Grazing, 'GrazingPerZ', GrazingPerZ)
        return GrazingPerZ



class ZooType:
    def __init__(self, stdpars, slicedpars):
        # zooplankton
        self.moZ = checkreplaceparam(stdpars, slicedpars, 'moZ')
        self.muZ = checkreplaceparam(stdpars, slicedpars, 'muZ')
        # grazing params
        self.Kp = checkreplaceparam(stdpars, slicedpars, 'Kp')
        self.deltaZ = checkreplaceparam(stdpars, slicedpars, 'deltaZ')

        self.muIntGraze = checkreplaceparam(stdpars, slicedpars, 'muIntGraze')
        self.kIntGraze = checkreplaceparam(stdpars, slicedpars, 'kIntGraze')

        self.pred = checkreplaceparam(stdpars, slicedpars, 'pred')
        self.deltaLambda = checkreplaceparam(stdpars, slicedpars, 'deltaLambda')

        self.pfn = stdpars['pfun_num']
        self.zn = stdpars['zoo_num']

        self.Vmax = self.muZ
        self.beta = 2
        self.ksat = self.Kp

        self.phylist = ['P' + str(i + 1) for i in range(self.pfn)]
        print(self.phylist)
        self.feedpref = [checkreplaceparam(stdpars, slicedpars, string) for string in self.phylist]

        print('feedpref',self.feedpref, self.muZ)

    def zoomortality(self, Z):
        # linear mortality of zooplankton, quadratic loss via higher order pred (below)
        total_moZ = self.moZ * Z
        return total_moZ

    def zoofeeding(self, P, func='anderson'):
        if func == 'anderson':
            Frho = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            GrazingProb = self.muZ / (self.ksat ** 2 + Frho)
            return GrazingProb

        elif func == 'fasham':
            Frho = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)]) # active switching coefficient
            Fp = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])  # total food available
            GrazingProb = self.muZ * (1 / (self.ksat * Fp + Frho))
            return GrazingProb

        elif func == 'vallina':
            Frho = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])  # active switching coefficient
            Fp = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])  # total food available
            GrazingProb = self.Vmax * (1 / Frho) * ((Fp ** self.beta) / (self.ksat ** self.beta + Fp ** self.beta))
            return GrazingProb

        else:
            print('no grazing formulation given, wrong func key')


    def fashamGP(self, P):
        Frho = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
        Fp = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])  # total food available

        GrazingProb = self.muZ * (1 / (self.ksat * Fp + Frho))
        # rhoP1 * P1 ** 2 + rhoP2 * P2 ** 2)
        # Gj = (self.muZ * rhoP1 * P1**2)/(kZ**2 + rhoP1 * P1**2 + rhoP2 * P2**2) * Zi
        return GrazingProb

    def andersonGP(self, P):
        # rhoP1 * P1 ** 2
        Frho = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
        GrazingProb = self.muZ / (self.ksat ** 2 + Frho)
        # rhoP1 * P1 ** 2 + rhoP2 * P2 ** 2)
        # Gj = (self.muZ * rhoP1 * P1**2)/(kZ**2 + rhoP1 * P1**2 + rhoP2 * P2**2) * Zi
        return GrazingProb

    def vallinaGP(self, P):
        Frho = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])  # active switching coefficient
        Fp = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])  # total food available
        # now calculate grazing probability of either Zoo-Type
        # print('p1', P[0], 'p2', P[1], 'Frho', Frho, 'Fp', Fp, 'Fq/Frho', Fq/Frho)

        Gj = self.Vmax * (1 / Frho) * ((Fp ** self.beta) / (self.ksat ** self.beta + Fp ** self.beta))

        return Gj

    def fullgrazing(self, GjJ, P, Zi):
        # phytouptake per zooplankton for each phyto
        Iprob = [GjJ * (self.feedpref[i] * P[i] ** 2) for i in range(self.pfn)]  # grazeprob per each PFT
        #perhaps here = multiply by Z before summing Iprobs!
        Itots = sum([Iprob[i] for i in range(self.pfn)])
        Itot = Itots * Zi
        return Itot

    def assimgrazing(self, Itoti):
        AssimGrazing = self.deltaZ * Itoti #* Zi
        return AssimGrazing

    def unassimilatedgrazing(self, Itoti):
        UnAsGraze = (1. - self.deltaZ) * Itoti #* Zi
        return UnAsGraze

    def interzoograze(self,i, Z):
        totalgraze = self.muIntGraze * Z[0] / (Z[0] + self.kIntGraze) * Z[1]
        if i == 0:
            return -totalgraze
        if i == 1:
            return self.deltaLambda * totalgraze
        else:
            return 0

    def unassiminterzoogr(self,i,Z):
        totalgraze = self.muIntGraze * Z[0] / (Z[0] + self.kIntGraze) * Z[1]
        if i == 0:
            return 0
        if i == 1:
            return (1 - self.deltaLambda) * totalgraze
        else:
            return 0

    def higherorderpred(self, Zi):
        total_moZ = self.pred * Zi ** 2  # sum(Z)
        return total_moZ / self.zn
