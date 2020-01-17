#!/usr/bin/env python
# -*- coding: utf-8 -*-

#import pandas
import numpy as np
#import scipy.interpolate as intrp
from PhytoMFTM.AuxFuncs import sliceparams, sliceoffparams, checkreplaceparam
from PhytoMFTM.ModelClasses import Forcing
from lmfit import minimize, Parameters#, Parameter, report_fit



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
        #tdp = np.exp(0.063 * Tmld)
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
#        self.interfeedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoointlistfeed]
#        self.intergrazedpref = [checkreplaceparam(allpars, slicedpars, string) for string in self.zoointlistgrazed]

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
    """"""
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

def cariaco(x,t,modelsetup, q):
    """System of ODEs"""

    N, P, Z, D, outputlist = modelsetup.timestep_init(x)
    #print(N,P,Z,D)
    #print(x)
    # N = [Ni]
    # P = [P1]
    # Z = [Z1]
    # D = [D]

    MLD = physx.MLD(t)  # MLD = [int_MLD, deriv_MLD]
    N0 = physx.N0(t)  # N0 = [Ni0,P0,Si0]
    #Si0 = physx.Si0(t)
    PAR = physx.PAR(t)
    Tmld = physx.Tmld(t)

    Mix = physx.omegaMix(MLD)  # i.e. there is constant mixing & increased mix when MLD shallowing
    #print('Mix',Mix)
    Mix_D = physx.omegaMix(MLD, type='D')  # i.e. there is constant mixing & increased mix when MLD shallowing

    # Grazing
    Gj = z.zoofeeding(P, Z, D, func='hollingtypeIII')  # feeding probability for all food
    ZooFeeding = z.fullgrazing(Gj, P, Z, D)

    PTempDepGrow = p.tempdepgrowth(Tmld)
    PNutUptake = p.uptake(N)
    PLightHarv = p.lightharvesting(MLD[0], PAR, P, sum(PTempDepGrow)) * 24/75  # (C to Chl)
    # Phytoplankton Fluxes
    PGains = PTempDepGrow * PNutUptake * PLightHarv * P

    PLinMort = p.mortality(P, type='linear')
    PQuadMort = p.mortality(P, type='quadratic')
    PMortality = PLinMort + PQuadMort
    PZooGrazed = p.zoograzing(Gj, P, Z, D)
    PMixing = P * Mix
    PLosses = PZooGrazed + PMortality + PMixing

    # Zooplankton Fluxes
    ZGains = z.assimgrazing(ZooFeeding)
    ZLinMort = z.mortality(Z, type='linear')
    ZQuadMort = z.mortality(Z, type='quadratic')
    ZMixing = Z * Mix
    ZLosses = ZLinMort + ZQuadMort + ZMixing

    # Detritus Fluxes
    ZUnassimFeedDetritus = z.unassimilatedgrazing(ZooFeeding, pool='D')
    DGains = sum(ZUnassimFeedDetritus) + sum(ZLinMort) + sum(PMortality)
    DRemin = d.remineralisation(D)
    DZooGrazed = d.zoograzing(Gj, D, Z)
    DMixing = D * Mix_D
    DLosses = DZooGrazed + DRemin + DMixing

    ZUnassimFeedNitrate = z.unassimilatedgrazing(ZooFeeding, pool='N')
    NMixing = Mix * (N0 - N)

    Px = PGains - PLosses
    Nx = - sum(PGains) + DRemin + sum(ZUnassimFeedNitrate) + NMixing# Nutrient draw down
    Zx = ZGains - ZLosses  # Zooplankton losses due to mortality and mixing
    Dx = DGains - DLosses   # Detritus

    out = [Nx, Px, Zx, Dx]

    outputlist[0] = PTempDepGrow
    outputlist[1] = PNutUptake
    outputlist[2] = PLightHarv
    outputlist[3] = PGains

    outputlist[4] = PLinMort
    outputlist[5] = PQuadMort
    outputlist[6] = PMortality
    outputlist[7] = PZooGrazed
    outputlist[22] = PMixing
    outputlist[8] = PLosses

    outputlist[9] = ZGains

    outputlist[10] = ZLinMort
    outputlist[11] = ZQuadMort
    outputlist[12] = ZMixing
    outputlist[13] = ZLosses

    outputlist[14] = ZUnassimFeedDetritus
    outputlist[15] = DGains

    outputlist[16] = DRemin
    outputlist[17] = DZooGrazed
    outputlist[18] = DMixing
    outputlist[19] = DLosses

    outputlist[20] = NMixing
    outputlist[21] = ZUnassimFeedNitrate

    return np.concatenate([out,outputlist], axis=None) #,


# MODEL ODE
def empower(x,t,modelsetup, q):
    """System of ODEs"""

    N, P, Z, D, outputlist = modelsetup.timestep_init(x)
    #print(N,P,Z,D)
    #print(x)
    # N = [Ni]
    # P = [P1]
    # Z = [Z1]
    # D = [D]

    MLD = physx.MLD(t)  # MLD = [int_MLD, deriv_MLD]
    N0 = physx.N0(t)  # N0 = [Ni0,P0,Si0]
    PAR = physx.PAR(t)
    Tmld = physx.Tmld(t)

    Mix = physx.omegaMix(MLD)  # i.e. there is constant mixing & increased mix when MLD shallowing
    #print('Mix',Mix)
    Mix_D = physx.omegaMix(MLD, type='D')  # i.e. there is constant mixing & increased mix when MLD shallowing

    # Grazing
    Gj = z.zoofeeding(P, Z, D, func='hollingtypeIII')  # feeding probability for all food
    ZooFeeding = z.fullgrazing(Gj, P, Z, D)

    PTempDepGrow = p.tempdepgrowth(Tmld)
    PNutUptake = p.uptake(N)
    PLightHarv = p.lightharvesting(MLD[0], PAR, P, sum(PTempDepGrow)) * 24/75  # (C to Chl)
    # Phytoplankton Fluxes
    PGains = PTempDepGrow * PNutUptake * PLightHarv * P

    PLinMort = p.mortality(P, type='linear')
    PQuadMort = p.mortality(P, type='quadratic')
    PMortality = PLinMort + PQuadMort
    PZooGrazed = p.zoograzing(Gj, P, Z, D)
    PMixing = P * Mix
    PLosses = PZooGrazed + PMortality + PMixing

    # Zooplankton Fluxes
    ZGains = z.assimgrazing(ZooFeeding)
    ZLinMort = z.mortality(Z, type='linear')
    ZQuadMort = z.mortality(Z, type='quadratic')
    ZMixing = Z * Mix
    ZLosses = ZLinMort + ZQuadMort + ZMixing

    # Detritus Fluxes
    ZUnassimFeedDetritus = z.unassimilatedgrazing(ZooFeeding, pool='D')
    DGains = sum(ZUnassimFeedDetritus) + sum(ZLinMort) + sum(PMortality)
    DRemin = d.remineralisation(D)
    DZooGrazed = d.zoograzing(Gj, D, Z)
    DMixing = D * Mix_D
    DLosses = DZooGrazed + DRemin + DMixing

    ZUnassimFeedNitrate = z.unassimilatedgrazing(ZooFeeding, pool='N')
    NMixing = Mix * (N0 - N)

    Px = PGains - PLosses
    Nx = - sum(PGains) + DRemin + sum(ZUnassimFeedNitrate) + NMixing# Nutrient draw down
    Zx = ZGains - ZLosses  # Zooplankton losses due to mortality and mixing
    Dx = DGains - DLosses   # Detritus

    out = [Nx, Px, Zx, Dx]

    outputlist[0] = PTempDepGrow
    outputlist[1] = PNutUptake
    outputlist[2] = PLightHarv
    outputlist[3] = PGains

    outputlist[4] = PLinMort
    outputlist[5] = PQuadMort
    outputlist[6] = PMortality
    outputlist[7] = PZooGrazed
    outputlist[22] = PMixing
    outputlist[8] = PLosses

    outputlist[9] = ZGains

    outputlist[10] = ZLinMort
    outputlist[11] = ZQuadMort
    outputlist[12] = ZMixing
    outputlist[13] = ZLosses

    outputlist[14] = ZUnassimFeedDetritus
    outputlist[15] = DGains

    outputlist[16] = DRemin
    outputlist[17] = DZooGrazed
    outputlist[18] = DMixing
    outputlist[19] = DLosses

    outputlist[20] = NMixing
    outputlist[21] = ZUnassimFeedNitrate

    return np.concatenate([out,outputlist], axis=None) #,



######### PARAMETER SETUP #############

# set up basic parameters of model:
standardparams = Parameters()

# number of phytoplankton func types
# standardparams.add('pfun_num', value=4, vary=False)
# number of zooplankton groups
# standardparams.add('zoo_num', value=2, vary=False)

# mld - related
standardparams.add('kappa', value=0.15, vary=False) # vary=False) # min=0.09, max=0.11) # Diffusive mixing across thermocline (m*d^-1)
standardparams.add('deltaD_N', value=0.2, vary=False)   # Nitrate Remineralization rate (d^-1)

standardparams.add('kw', value=0.04, vary=False)     # Light attenuation constant of water (m^-1)

standardparams.add('kc', value=0.03, vary=False)      # Light attenuation via phytoplankton pigment (m^-1)
standardparams.add('alpha', value=0.15, vary=False)  # initial slope of the P-I curve
standardparams.add('VpMax', value=1., vary=False)    # maximum photosynthetic rate

standardparams.add('v', value=0., vary=False)      # Sinking of Phytoplankton from Mixed Layer
standardparams.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)

# p - related
standardparams.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)

standardparams.add('ratioSi', value=0, vary=False)  # Silicate ratio ## FALLBACK PARAM FOR OTHER PFTs
standardparams.add('U_N', value=0, vary=False)    # Nitrate Half Saturation Constant
standardparams.add('U_Si', value=0, vary=False)   # Silicate Half Saturation Constant
standardparams.add('muP', value=0, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 1 (e.g. DIATOMS)
ptype1 = Parameters()
ptype1.add('pt1_ratioSi', value=1., vary=False)  # Silicate ratio
ptype1.add('pt1_U_Si', value=0.5, vary=False)   # Silicate Half Saturation Constant
ptype1.add('pt1_U_N', value=0.9, vary=False)    # Nitrate Half Saturation Constant
ptype1.add('pt1_muP', value=1.074, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 2 (e.g. Haptos)
ptype2 = Parameters()
#ptype2.add('pt2_OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)
ptype2.add('pt2_U_N', value=1.5, vary=False)    # Nitrate Half Saturation Constant
ptype2.add('pt2_muP', value=0.51, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 3 (e.g. Cyanos)
ptype3 = Parameters()
ptype3.add('pt3_U_N', value=2.5, vary=False)    # Nitrate Half Saturation Constant
ptype3.add('pt3_muP', value=.4, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 4 (e.g. Dinos)
ptype4 = Parameters()
ptype4.add('pt4_U_N', value=1.5, vary=False)    # Nitrate Half Saturation Constant
ptype4.add('pt4_muP', value=0.3, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 5 (e.g. Others)
ptype5 = Parameters()
ptype5.add('pt5_U_N', value=1.84, vary=False)    # Nitrate Half Saturation Constant
ptype5.add('pt5_muP', value=0.51, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# z - related
#z grazing related
standardparams.add('moZ', value=0.07, vary=False)        # Zooplankton mortality (d^-1)
standardparams.add('deltaZ', value=0.75, vary=False)    # Zooplankton Grazing assimilation coefficient (-)
standardparams.add('deltaLambda', value=0.75, vary=False)    # Zooplankton Inter-Grazing assimilation coefficient (-)
standardparams.add('muIntGraze', value=0.1, vary=False)  # InterZooGrazing maximum grazing rate
standardparams.add('kIntGraze', value=0.5, vary=False)  # InterZooGrazing saturation constant

standardparams.add('Kp', value=0, vary=False)     # Zooplankton Grazing saturation constant (-)
standardparams.add('pred', value=0, vary=False)  # quadratic higher order predation rate on zooplankton
standardparams.add('muZ', value=0, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# set up zooplankton type 1 (e.g. MIKRO zooplankton)
ztype1 = Parameters()
ztype1.add('zt1_muZ', value=0.5, min=.2, max=1.5)    # Zooplankton maximum grazing rate (d^-1)

ztype1.add('zt1_Kp', value=1., min=.2, max=1.5)       # Zooplankton Grazing saturation constant (-)
ztype1.add('zt1_pred', value=0.06, min=.001, max=0.2)    # quadratic higher order predation rate on zooplankton

# set up zooplankton type 2 (e.g. MESO zooplankton)
ztype2 = Parameters()
ztype2.add('zt2_muZ', value=0.5, min=.2, max=1.5)    # Zooplankton maximum grazing rate (d^-1)

ztype2.add('zt2_Kp', value=1.3, min=.2, max=1.5)       # Zooplankton Grazing saturation constant (-)
ztype2.add('zt2_pred', value=ztype1['zt1_pred'].value, vary=False)    # quadratic higher order predation rate on zooplankton

"""
add feeding pref params!
Z1P1
Z1P2
Z2P1
Z2P2
"""
# MIKRO
ztype1.add('zt1_P1', value=0, vary=False)  # Diatoms
ztype1.add('zt1_P2', value=0.25, vary=False)  # Hapto
ztype1.add('zt1_P3', value=0.25, vary=False)  # Cyano
ztype1.add('zt1_P4', value=0.25, vary=False)  # Dino
ztype1.add('zt1_P5', value=0.25, vary=False)  # Others
# MESO
ztype2.add('zt2_P1', value=1/6, vary=False)
ztype2.add('zt2_P2', value=1/6, vary=False)
ztype2.add('zt2_P3', value=1/6, vary=False)
ztype2.add('zt2_P4', value=1/6, vary=False)
ztype2.add('zt2_P5', value=1/6, vary=False)

# inter zoo feeding
# MIKRO
ztype1.add('zt1_Zint_feed1', value=0, vary=False)
ztype1.add('zt1_Zint_feed2', value=0, vary=False)
# MESO
ztype2.add('zt2_Zint_feed1', value=1/6, vary=False)
ztype2.add('zt2_Zint_feed2', value=0, vary=False)

# CONVERT FEEDPREFS TO GRAZEPREF FOR CALCULATION OF GRAZING
ztype1.add('zt1_Zint_grazed1', value=ztype1['zt1_Zint_feed1'].value, vary=False)
ztype1.add('zt1_Zint_grazed2', value=ztype2['zt2_Zint_feed1'].value, vary=False)

ztype2.add('zt2_Zint_grazed1', value=ztype1['zt1_Zint_feed2'].value, vary=False)
ztype2.add('zt2_Zint_grazed2', value=ztype2['zt2_Zint_feed2'].value, vary=False)

ptype1.add('pt1_Z1', value=ztype1['zt1_P1'].value, vary=False)
ptype1.add('pt1_Z2', value=ztype2['zt2_P1'].value, vary=False)

ptype2.add('pt2_Z1', value=ztype1['zt1_P2'].value, vary=False)
ptype2.add('pt2_Z2', value=ztype2['zt2_P2'].value, vary=False)

ptype3.add('pt3_Z1', value=ztype1['zt1_P3'].value, vary=False)
ptype3.add('pt3_Z2', value=ztype2['zt2_P3'].value, vary=False)

ptype4.add('pt4_Z1', value=ztype1['zt1_P4'].value, vary=False)
ptype4.add('pt4_Z2', value=ztype2['zt2_P4'].value, vary=False)

ptype5.add('pt5_Z1', value=ztype1['zt1_P5'].value, vary=False)
ptype5.add('pt5_Z2', value=ztype2['zt2_P5'].value, vary=False)

parameters = standardparams + ptype1 + ztype1 #+ ztype2 + ptype2 + ptype3 + ptype4 + ptype5

######### MODEL EVALUATION CODE #############

ms = ModelSetup(parameters)

n,p,z,d = ms.classes

physx = ms.physics


N0 = 5
P0 = 0.1
Z0 = 0.1
D0 = 0.1

initnut = [N0 for i in range(n.num)]
initphy = [P0 for i in range(p.num)]
initzoo = [Z0 for i in range(z.num)]
initdet = [D0 for i in range(d.num)]
initout = [0 for i in range(23)]
initcond = np.concatenate([initnut, initphy, initzoo, initdet,initout], axis=None)

timedays = np.arange(0., 5 * 365., 1.0)

import time
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# INTEGRATE:
tos = time.time()
print('starting integration')
outarray = odeint(empower, initcond, timedays, args=(ms,1))
tos1 = time.time()
print('finished after %4.3f sec' % (tos1 - tos))

#print(outarray)


numcols = 2
f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, numcols, sharex='col')#, sharey='row')


plt.setp((ax1, ax2, ax3, ax4), xticks=[1,60,120,180,240,300,365])
from matplotlib.ticker import MaxNLocator
for axe in (ax1, ax2, ax3, ax4):
    for i in range(numcols):
        axe[i].get_yaxis().set_major_locator(MaxNLocator(nbins=4))
        axe[i].tick_params(top=True, right=True)

# PLOTTING
timedays_ly = timedays[1:366]
outarray_ly = outarray[1460:1825]

# color vectors
#colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', 'grey']
alphas = [1., 0.8, 0.6, 0.4]
lws = [2, 2.5, 4, 5.5]

ax1[0].set_title('model output')

dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
dpm_cumsum = np.cumsum(dayspermonth) - np.array(dayspermonth)/2 #- 15
print(timedays_ly)


# Figure 1
# N
N_Max = np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) + np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) * 0.1
ax1[0].scatter(dpm_cumsum, ms.physics.forcing.verif.N, label='WOA data')
ax1[0].plot(timedays_ly, outarray_ly[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
ax1[0].set_ylabel('Nutrients \n' '[µM N]', multialignment='center', fontsize=10)
ax1[0].set_ylim(0, N_Max)
ax1[0].legend(fontsize='x-small')


ChlConv = True
# Phyto
CtoChla = 75  # g/g
MolarMassC = 12.0107
CtoNratioPhyto = 6.625
muMolartoChlaconvfactor = CtoChla / MolarMassC / CtoNratioPhyto  # Chla as mg/m-3 to

ax2[0].scatter(dpm_cumsum, np.array(ms.physics.forcing.verif.chla) * muMolartoChlaconvfactor, label='MODISaq data')


Pall = outarray_ly[:,1]
P_Max = np.max(Pall) + 0.9 * np.max(Pall)

ax2[0].plot(timedays_ly, Pall, c=colors[4], lw=lws[1], label='Model')
ax2[0].legend(fontsize='x-small')
ax2[0].set_ylabel('Phytoplankton \n' '[µM N]', multialignment='center', fontsize=10)
ax2[0].set_ylim(0, P_Max)

# Z
Zall = outarray_ly[:,2]
Z_Max = np.max(Zall) + 0.1 * np.max(Zall)

ax3[0].plot(timedays_ly, Zall, c=colors[4], lw=lws[1])
ax3[0].set_ylabel('Zooplankton \n' '[µM N]', multialignment='center', fontsize=9)
ax3[0].tick_params('y', labelsize=10)
ax3[0].set_ylim(0, Z_Max)
#ax4[i_plot].set_title('Zooplankton')

D_Max = np.max(outarray_ly[:, 3]) + 0.2 * np.max(outarray_ly[:, 3])
# D
ax4[0].plot(timedays_ly, outarray_ly[:, 3], c=colors[1], lw=lws[0], alpha=alphas[0])
ax4[0].set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)
ax4[0].set_ylim(0,D_Max)
ax4[0].set_xlabel('Day in year')
# Legend

## PHYSICS ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
muplot = 1
ax1[muplot].set_title('model forcing')

ax4[muplot].set_xlabel('Day in year')

NOX = ms.physics.forcing.NOX.return_interpvalattime(timedays_ly)
NOXdat = ms.physics.forcing.NOX.forcingfile
print(NOX)
print(NOXdat)
ax1[muplot].plot(timedays_ly, NOX, c=colors[5], lw=lws[0], alpha=alphas[0])
ax1[muplot].scatter(dpm_cumsum, NOXdat[0:12], c=colors[5])
ax1[muplot].set_ylabel('$N_0$ \n' '[µM]', multialignment='center', fontsize=10)
ax1[muplot].set_ylim(0., N_Max)
#ax1[muplot].invert_yaxis()

MLD = ms.physics.forcing.MLD.return_interpvalattime(timedays_ly)
MLDdat = ms.physics.forcing.MLD.forcingfile
MLD_max = np.max(MLD) + 0.1 * np.max(MLD)
ax2[muplot].plot(timedays_ly, MLD, c=colors[5], lw=lws[0], alpha=alphas[0])
ax2[muplot].scatter(dpm_cumsum, MLDdat[0:12], c=colors[5])
ax2[muplot].set_ylabel('MLD \n' '[m]', multialignment='center', fontsize=10)
ax2[muplot].set_ylim(0, MLD_max) # 400 for biotrans, 100 for Papa
ax2[muplot].invert_yaxis()

PAR = ms.physics.forcing.PAR.return_interpvalattime(timedays_ly)
PARdat = ms.physics.forcing.PAR.forcingfile
PAR_max = np.max(PAR) + 0.1 * np.max(PAR)
ax3[muplot].plot(timedays_ly, PAR, c=colors[5], lw=lws[0], alpha=alphas[0])
ax3[muplot].scatter(dpm_cumsum, PARdat[0:12], c=colors[5])
ax3[muplot].set_ylabel('PAR \n' '[E $m^{−2}$ $s^{−1}$]', multialignment='center', fontsize=10)
ax3[muplot].set_ylim(0, PAR_max)
# ax1[muplot].invert_yaxis()

Tmld = ms.physics.forcing.SST.return_interpvalattime(timedays_ly)
Tmlddat = ms.physics.forcing.SST.forcingfile
Tmld_max = np.max(Tmld) + 0.1 * np.max(Tmld)
ax4[muplot].plot(timedays_ly, Tmld, c=colors[5], lw=lws[0], alpha=alphas[0])
ax4[muplot].scatter(dpm_cumsum, Tmlddat[0:12], c=colors[5])
ax4[muplot].set_ylabel('$T_{MLD}$ \n' '[°C]', multialignment='center', fontsize=10)
ax4[muplot].set_ylim(0, Tmld_max)
# ax1[muplot].invert_yaxis()

# Defining custom 'xlim' and 'ylim' values.
xlim = (0, 365)

# Setting the values for all axes.
plt.setp((ax1, ax2, ax3, ax4), xlim=xlim)

f1.align_ylabels()

plt.subplots_adjust(hspace=0.1)

plt.tight_layout()
plt.show()