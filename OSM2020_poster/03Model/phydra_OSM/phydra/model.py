#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

# TODO:
#  - import modelsetup (not necessary?)

def cariaco(x,t,modelsetup, q):
    """System of ODEs"""
    # TODO: add sinking parameter to all components, instead of only for Detritus? or keep it like dat?

    N, P, Z, D, outputlist = modelsetup.timestep_init(x)

    n, p, z, d = modelsetup.classes

    physx = modelsetup.physics

    #print(N,P,Z,D)
    #print(x)
    # N = [Ni]
    # P = [P1]
    # Z = [Z1]
    # D = [D]

    X258 = physx.X258(t)  # X258 = [int_MLD, deriv_MLD]
    N0 = physx.N0(t)  # N0 = [Ni0,P0,Si0]
    #Si0 = physx.Si0(t)
    PAR = physx.PAR(t)
    Tmld = physx.Tmld(t)

    Mix = physx.wMixMix(X258)  # i.e. there is constant mixing & increased mix when MLD shallowing
    #print('Mix',Mix)
    Mix_D = physx.wMixMix(X258, type='D')  # i.e. there is constant mixing & increased mix when MLD shallowing

    # Grazing
    Gj = z.zoofeeding(P, Z, D, func='hollingtypeIII')  # feeding probability for all food
    ZooFeeding = z.fullgrazing(Gj, P, Z, D)

    # TODO: make sure tempdepgrowth is correct implementation, not empower style
    PTempDepGrow = p.tempdepgrowth(Tmld)
    PNutUptake = p.uptake(N)
    # TODO: Update the light limiting term to use Steele instead of Smith, also depth of box is 100m not x258
    PLightHarv = p.lightharvesting(X258[0], PAR, P, sum(PTempDepGrow)) * 24/75  # (C to Chl)
    # Phytoplankton Fluxes
    PGains = PTempDepGrow * PNutUptake * PLightHarv * P

    PLinMort = p.mortality(P, type='linear')
    # TODO: remove PQuadMort
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
    n, p, z, d = modelsetup.classes
    physx = modelsetup.physics

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

    return np.concatenate([out,outputlist], axis=None)