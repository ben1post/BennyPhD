#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def cariaco(x,t,modelsetup, q):
    """System of ODEs - model structure adapted to CARIACO time series setting"""
    N, P, Z, D, outputlist = modelsetup.timestep_init(x)

    n, p, z, d = modelsetup.classes

    physx = modelsetup.physics

    #print(N,P,Z,D)
    #print(x)
    #print(type(P))
    # N = [Ni]
    # P = [P1]
    # Z = [Z1]
    # D = [D]

    X258 = physx.X258(t)  # X258 = [int_X258, deriv_X258]
    Ni0 = physx.N0(t) # N0 = [Ni0,P0,Si0]
    Si0 = physx.Si0(t)
    NUT0 = [Ni0, Si0]
    PAR = physx.PAR(t)
    Tmld = physx.Tmld(t)

    Mix = physx.wMix(X258)  # i.e. there is constant mixing & increased mix when MLD shallowing
    #print('Mix',Mix)
    if X258[0] < 100: Mix = Mix * 1.5

    if X258[0] < 60: Mix = Mix * 2.5

    # Grazing
    Gj = z.zoofeeding(P, Z, D, func='hollingtypeIII')  # feeding probability for all food
    ZooFeeding = z.fullgrazing(Gj, P, Z, D)

    PGrowth = p.growth()
    PTempDepGrow = p.tempdepgrowth(Tmld)
    PNutUptake = p.uptake(N)
    PLightHarv = p.lightharvesting(physx.BoxDepth, PAR, P)
    # Phytoplankton Fluxes
    PGains = PGrowth * PTempDepGrow * PNutUptake * PLightHarv * P

    PMortality = p.mortality(P, type='linear')
    PZooGrazed = p.zoograzing(Gj, P, Z, D)
    PSinking = p.sinking(physx.BoxDepth, P)
    PLosses = PZooGrazed + PMortality + PSinking

    # Zooplankton Fluxes
    ZGains = z.assimgrazing(ZooFeeding)
    ZLinMort = z.mortality(Z, type='linear')
    ZQuadMort = z.mortality(Z, type='quadratic')

    #if t > 365*3 + 365*5:
    #    ZQuadMort = ZQuadMort/2

    ZLosses = ZLinMort + ZQuadMort

    # Detritus Fluxes
    ZUnassimFeedDetritus = z.unassimilatedgrazing(ZooFeeding, pool='D')
    DGains = sum(ZUnassimFeedDetritus) + sum(ZLinMort) + sum(PMortality)
    DRemin = d.remineralisation(D)
    DZooGrazed = d.zoograzing(Gj, D, Z)
    DSinking = d.sinking(D, physx.BoxDepth)
    DLosses = DZooGrazed + DRemin + DSinking

    ZUnassimFeedNitrate = z.unassimilatedgrazing(ZooFeeding, pool='N')
    NUTMixing = n.mixing(NUT0, N, Mix)  # Mix * (N0 - N)
    NMixing = NUTMixing[0]
    SiMixing = NUTMixing[1]

    Px = PGains - PLosses
    Nix = - sum(PGains) + DRemin + sum(ZUnassimFeedNitrate) + NMixing# Nutrient draw down
    Six = - PGains[0] + SiMixing
    Nx = [Nix, Six]
    Zx = ZGains - ZLosses  # Zooplankton losses due to mortality and mixing
    Dx = DGains - DLosses   # Detritus

    out = np.concatenate((Nx, Px, Zx, Dx), axis=None)

    outputlist[0] = sum(PTempDepGrow)
    outputlist[1] = sum(PNutUptake)
    outputlist[2] = sum(PLightHarv)
    outputlist[3] = sum(PGains)

    outputlist[4] = Mix #0 #PLinMort
    outputlist[5] = 0 #PQuadMort
    outputlist[6] = sum(PMortality)
    outputlist[7] = sum(PZooGrazed)
    outputlist[22] = sum(PSinking) #PMixing
    outputlist[8] = sum(PLosses)

    outputlist[9] = ZGains

    outputlist[10] = ZLinMort
    outputlist[11] = ZQuadMort
    outputlist[12] = SiMixing[0] # 0 # ZMixing
    outputlist[13] = ZLosses

    outputlist[14] = ZUnassimFeedDetritus
    outputlist[15] = DGains

    outputlist[16] = DRemin
    outputlist[17] = DZooGrazed
    outputlist[18] = DSinking #DMixing
    outputlist[19] = DLosses

    outputlist[20] = NMixing[0]
    outputlist[21] = ZUnassimFeedNitrate

    return np.concatenate((out, outputlist), axis=None)


def cariacoNPZD(x,t,modelsetup, q):
    """System of ODEs - model structure adapted to CARIACO time series setting"""
    N, P, Z, D, outputlist = modelsetup.timestep_init(x)

    n, p, z, d = modelsetup.classes

    physx = modelsetup.physics

    #print(N,P,Z,D)
    #print(x)
    #print(type(P))
    # N = [Ni]
    # P = [P1]
    # Z = [Z1]
    # D = [D]

    X258 = physx.X258(t)  # X258 = [int_X258, deriv_X258]
    Ni0 = physx.N0(t) # N0 = [Ni0,P0,Si0]
    #Si0 = physx.Si0(t)
    NUT0 = Ni0#, Si0]
    PAR = physx.PAR(t)
    Tmld = physx.Tmld(t)

    Mix = physx.wMix(X258)  # i.e. there is constant mixing & increased mix when MLD shallowing
    #print('Mix',Mix)
    #if X258[0] < 100: Mix = Mix * 1.5

    #if X258[0] < 60: Mix = Mix * 2.5

    # Grazing
    Gj = z.zoofeeding(P, Z, D, func='hollingtypeIII')  # feeding probability for all food
    ZooFeeding = z.fullgrazing(Gj, P, Z, D)

    PGrowth = p.growth()
    PTempDepGrow = p.tempdepgrowth(Tmld)
    PNutUptake = p.uptake(N)
    PLightHarv = p.lightharvesting(physx.BoxDepth, PAR, P)
    # Phytoplankton Fluxes
    PGains = PGrowth * PTempDepGrow * PNutUptake * PLightHarv * P

    PMortality = p.mortality(P, type='linear')
    PZooGrazed = p.zoograzing(Gj, P, Z, D)
    PSinking = p.sinking(physx.BoxDepth, P)
    PLosses = PZooGrazed + PMortality + PSinking

    # Zooplankton Fluxes
    ZGains = z.assimgrazing(ZooFeeding)
    ZLinMort = z.mortality(Z, type='linear')
    ZQuadMort = z.mortality(Z, type='quadratic')

    #if t > 365*3 + 365*5:
    #    ZQuadMort = ZQuadMort/2

    ZLosses = ZLinMort + ZQuadMort

    # Detritus Fluxes
    ZUnassimFeedDetritus = z.unassimilatedgrazing(ZooFeeding, pool='D')
    DGains = sum(ZUnassimFeedDetritus) + sum(ZLinMort) + sum(PMortality)
    DRemin = d.remineralisation(D)
    DZooGrazed = d.zoograzing(Gj, D, Z)
    DSinking = d.sinking(D, physx.BoxDepth)
    DLosses = DZooGrazed + DRemin + DSinking

    ZUnassimFeedNitrate = z.unassimilatedgrazing(ZooFeeding, pool='N')
    NUTMixing = n.mixing(NUT0, N, Mix)  # Mix * (N0 - N)
    NMixing = NUTMixing#[0]
   # SiMixing = NUTMixing[1]

    Px = PGains - PLosses
    Nix = - sum(PGains) + DRemin + sum(ZUnassimFeedNitrate) + NMixing# Nutrient draw down
    #Six = - PGains[0] + SiMixing
    Nx = Nix#, Six]
    Zx = ZGains - ZLosses  # Zooplankton losses due to mortality and mixing
    Dx = DGains - DLosses   # Detritus

    out = np.concatenate((Nx, Px, Zx, Dx), axis=None)

    outputlist[0] = sum(PTempDepGrow)
    outputlist[1] = sum(PNutUptake)
    outputlist[2] = sum(PLightHarv)
    outputlist[3] = sum(PGains)

    outputlist[4] = Mix #0 #PLinMort
    outputlist[5] = 0 #PQuadMort
    outputlist[6] = sum(PMortality)
    outputlist[7] = sum(PZooGrazed)
    outputlist[22] = sum(PSinking) #PMixing
    outputlist[8] = sum(PLosses)

    outputlist[9] = ZGains

    outputlist[10] = ZLinMort
    outputlist[11] = ZQuadMort
    outputlist[12] = 0#SiMixing[0] # 0 # ZMixing
    outputlist[13] = ZLosses

    outputlist[14] = ZUnassimFeedDetritus
    outputlist[15] = DGains

    outputlist[16] = DRemin
    outputlist[17] = DZooGrazed
    outputlist[18] = DSinking #DMixing
    outputlist[19] = DLosses

    outputlist[20] = NMixing[0]
    outputlist[21] = ZUnassimFeedNitrate

    return np.concatenate((out, outputlist), axis=None)


# MODEL ODE
def empower(x,t,modelsetup, q):
    """System of ODEs - model setup for EMPOWER-like NPZD model implementation"""

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