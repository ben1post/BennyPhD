#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from phydra.aux import checkreplaceparam

# TODO:
#  - make sure all CARIACO functions are represented here
#  -

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
        self.moP_quad = checkreplaceparam(allpars, slicedpars, 'moP_quad')

        self.ratioSi = checkreplaceparam(allpars, slicedpars, 'ratioSi')

        self.pfn = allpars['phyto_num'].value
        self.zn = allpars['zoo_num'].value

        self.kc = checkreplaceparam(allpars, slicedpars, 'kc')
        self.alpha = checkreplaceparam(allpars, slicedpars, 'alpha')
        self.VpMax = checkreplaceparam(allpars, slicedpars, 'VpMax')

        self.zoolist = ['Z'+str(j+1) for j in range(self.zn)]
        self.grazepref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoolist]

    def growth(self):
        return self.muP

    def uptake(self, N):
        """Michealis-Menten uptake of Nutrients"""
        # TODO: add functionality of multiple Nutrients/formulations
        Uptake = N[0] / (N[0] + self.U_N)
        return Uptake

    def lightharvesting(self, MLD, PAR, P, VPT, type='Smith'):
        """Light - modification of phytoplankton growth rate"""
        if type == 'Steele':
            lighthrv = 1. / (self.kw * MLD) * \
                       (-np.exp(1. - PAR / self.OptI) - (
                           -np.exp((1. - (PAR * np.exp(-self.kw * MLD)) / self.OptI))))
            return lighthrv

        # TODO: convert P to chlorophyll!
        if type == 'Smith':
            PAR = PAR
            kPAR = self.kw + self.kc * sum(P)
            x_0 = self.alpha * PAR * np.exp(- kPAR * 0)
            x_H = self.alpha * PAR * np.exp(- kPAR * MLD)
            VpH = (VPT / (kPAR * MLD)) * \
                  np.log(
                      (x_0 + np.sqrt(VPT ** 2 + x_0 ** 2)) /
                      (x_H + np.sqrt(VPT ** 2 + x_H ** 2))
                  )
            return VpH/VPT

    def tempdepgrowth(self, Tmld):
        """"""
        # tdp = np.exp(0.063 * Tmld)
        tdp = self.VpMax * 1.066 ** Tmld
        return tdp

    def mortality(self, P, type='linear'):
        """"""
        if type == 'linear':
            mortal = self.moP * P[self.num]
            return mortal
        elif type == 'quadratic':
            mortal = self.moP_quad * P[self.num] ** 2
            return mortal

    def zoograzing(self, Gj, P, Z, D):
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
        self.dn = allpars['det_num'].value

        self.beta = 2  # for Vallina KTW Grazing, beta = 1 : Holling Type II, beta = 2 : Holling Type III
        self.ksat = self.Kp

        self.phylist = ['P'+ str(i + 1) for i in range(self.pfn)] # + str(i + 1)
        self.detlist = ['D'+ str(i + 1) for i in range(self.dn)] # + str(i + 1)
        self.zoointlistfeed = ['Zint_feed' + str(j + 1) for j in range(self.zn)]  # list of feeding
        self.zoointlistgrazed = ['Zint_grazed' + str(j + 1) for j in range(self.zn)]
        # print(self.phylist, self.zoointlistfeed, self.zoointlistgrazed)

        self.feedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.phylist]
        self.detfeedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.detlist]
        print('D feedpref',self.detfeedpref,'P feedpref',self.feedpref)
        # self.interfeedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoointlistfeed]
        # self.intergrazedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoointlistgrazed]

        self.beta_feed = checkreplaceparam(allpars, slicedpars, 'beta_feed')
        self.kN_feed = checkreplaceparam(allpars, slicedpars, 'kN_feed')

    def zoofeeding(self, P, Z, D, func='anderson'):
        if func == 'anderson':
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            #FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP #+ FrhoZ
            GrazingProb = self.muZ / (self.ksat ** 2 + Frho)
            return GrazingProb

        elif func == 'fasham': #put the holling type num here
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            #FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP #+ FrhoZ  # active switching coefficient
            FpP = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])
            #FpZ = sum([self.interfeedpref[j] * Z[j] for j in range(self.zn)])
            Fp = FpP #+ FpZ  # total food available
            GrazingProb = self.muZ * (1 / (self.ksat * Fp + Frho))
            return GrazingProb

        elif func == 'vallina':
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            #FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP #+ FrhoZ  # active switching coefficient
            FpP = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])
            #FpZ = sum([self.interfeedpref[j] * Z[j] for j in range(self.zn)])
            Fp = FpP #+ FpZ  # total food available
            GrazingProb = self.muZ * (1 / Frho) * ((Fp ** self.beta) / (self.ksat ** self.beta + Fp ** self.beta))
            return GrazingProb

        elif func == 'hollingtypeIII':
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            FrhoD = sum([self.detfeedpref[j] * D[j] ** 2 for j in range(self.dn)])
            Frho = FrhoP + FrhoD
            #print('Frho',Frho,'FrhoP',FrhoP,P,'FrhoD',FrhoD,D)
            GrazingProb = self.muZ / (self.ksat ** 2 + Frho)
            return GrazingProb

        else:
            print('no grazing formulation given, wrong func key')

    def fullgrazing(self, Gj, P, Z, D):
        # phytouptake + zooplankton per zooplankton for each phyto
        IprobP = [Gj[self.num] * (self.feedpref[i] * P[i] ** 2) for i in range(self.pfn)]  # grazeprob per each PFT
        IprobD = [Gj[self.num] * (self.detfeedpref[j] * D[j] ** 2) for j in range(self.dn)]
        Iprob = IprobP + IprobD
        #perhaps here = multiply by Z before summing Iprobs!
        Itots = sum(Iprob)
        Itot = Itots * Z[self.num]
        #print('Itots', Itots,Z, 'IprobP', IprobP, P, 'IprobD', IprobD, D)
        return Itot

    def assimgrazing(self, ZooFeeding):
        # AssimGrazing = self.deltaZ * ZooFeeding[self.num]
        AssimGrazing = self.beta_feed * self.kN_feed * ZooFeeding[self.num]
        return AssimGrazing

    def unassimilatedgrazing(self, ZooFeeding, pool='N'):
        #UnAsGraze = (1. - self.deltaZ) * ZooFeeding[self.num]
        if pool == 'N':
            UnAsGraze = self.beta_feed * (1-self.kN_feed) * ZooFeeding[self.num]
            return UnAsGraze
        elif pool == 'D':
            UnAsGraze = (1-self.beta_feed) * ZooFeeding[self.num]
            return UnAsGraze

    def mortality(self, Z, type='linear'):
        if type == 'linear':
            total_moZ = self.moZ * Z[self.num]
            return total_moZ
        if type == 'quadratic':
            total_moZ = self.pred * Z[self.num] ** 2
            return total_moZ


class Detritus:
    def __init__(self, allpars, slicedpars, num):
        self.num = num
        self.deltaD = checkreplaceparam(allpars, slicedpars, 'deltaD_N')

        self.zn = allpars['zoo_num'].value
        self.zoolist = ['Z' + str(j + 1) for j in range(self.zn)]
        self.detgrazepref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoolist]

    def remineralisation(self, D):
        return D[self.num] * self.deltaD

    def zoograzing(self, Gj, D, Z):
        """"""
        # take the general grazing term from each zooplankton, multiply by phyto fraction and sum
        Grazing = [Gj[j] * (self.detgrazepref[j] * D[self.num] ** 2) * Z[j] for j in range(self.zn)]
        GrazingPerZ = sum(Grazing)
        return GrazingPerZ