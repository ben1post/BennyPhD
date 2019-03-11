#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import pandas
from PhytoMFTM.AuxFuncs import zoograzeloss_phy,zoogrowth_zoo,grazing, dailyinterp, firstderivspl



# parameters for interpolation
kmld = 3
smld = 0
kindmld = "spline"
kn0x = 5
sn0x = None
kindn0x = "spline"
kpar = 5
spar = None
kindpar = "spline"
ksst = 5
ssst = None
kindsst = "spline"
ksi0x = 5
ssi0x = None
kindsi0x = "spline"

# read environmental forcings
MLDfile = pandas.read_csv('Forcing/MLD')
MLD = list(MLDfile['press'])
MLD.append(MLDfile['press'][0])

NOXfile = pandas.read_csv('Forcing/N0')
NOX = list(NOXfile['as.numeric(NO3_NO2)'])
NOX.append(NOXfile['as.numeric(NO3_NO2)'][0])

SiOXfile = pandas.read_csv('Forcing/Si0')
SiOX = list(SiOXfile['as.numeric(Silicate)'])
SiOX.append(SiOXfile['as.numeric(Silicate)'][0])

SSTfile = pandas.read_csv('Forcing/SST')
SST = list(SSTfile['temp'])
SST.append(SSTfile['temp'][0])

PARfile = pandas.read_csv('Forcing/PAR')
PAR = list(PARfile['value'])
PAR.append(PARfile['value'][0])



def simpleN2P2ZD(x, t, paras, pClass, zClass):
    N = x[0] # Nitrate
    Si = x[1] # Silicate
    D = x[2] # Detritus

    zn = paras['zoo_num'].value
    Z = [x[3 +i] for i in range(zn)] # Zooplankton

    z = zClass

    pfn = paras['pfun_num'].value
    P = [x[3 +zn +i] for i in range(pfn)] # Phytoplankton

    p = pClass
    # Interpolations of Forcings
    int_MLD = dailyinterp(MLD, t, kind=kindmld, k=kmld, s=smld)
    int_NOX = dailyinterp(NOX, t, kind=kindn0x, k=kn0x, s=sn0x)
    if int_NOX < 0 : int_NOX = 0  # do not allow negative Nitrate values
    int_SIOX = dailyinterp(SiOX, t, kind=kindsi0x, k=ksi0x, s=ssi0x)
    int_PAR = dailyinterp(PAR, t, kind=kindpar, k=kpar, s=spar)
    int_SST = dailyinterp(SST, t, kind=kindsst, k=ksst, s=ssst)

    # Derivatives of Forcings
    deriv_MLD = firstderivspl(MLD, t, k=kmld, s=smld)

    ## Non-Phytoplankton related processes
    # Mixing Processes
    DiffusiveMixing = (paras['kappa'].value + max(deriv_MLD, 0)) / int_MLD # i.e. there is constant mixing & increased loss with MLD shallowing
    ActiveMixing = deriv_MLD / int_MLD # i.e. concentration varies with the MLD depth

    # Remineralisation
    NRemineralization = paras['deltaD_N'].value * D
    SiRemineralization = paras['deltaD_Si'].value * D

    # Detritus
    DetritusMixing = D * DiffusiveMixing

    # Nutrient Mixing
    NMixing = DiffusiveMixing * (int_NOX - N)
    SiMixing = DiffusiveMixing * (int_SIOX - Si)

    # Zooplankton
    ZooMortality = [z[i].zoomortality(Z[i]) for i in range(zn)]

    y0 = NRemineralization + NMixing   # Nitrate upwelling and remineralisation
    y1 = SiRemineralization + SiMixing # Silicate upwelling and remineralisation
    y2 = sum(ZooMortality) - NRemineralization - SiRemineralization - DetritusMixing



    ## Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]

    # Light and Temperature
    LightHarvesting = [p[i].lightharvesting(int_MLD ,int_PAR) for i in range(pfn)]
    TemperatureDepGrowth = np.exp(0.063 * int_SST)

    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth)
             for i in range(pfn)]

    Ptot = sum(P)
    # Phytoplankton Losses
    GrazingMatrix = grazing(Ptot,Z, pfn,zn, z) # return all zooplankton feeding on all phytoplankton types

    PhytoGrazed = zoograzeloss_phy(GrazingMatrix, pfn ,zn) # sums up z.grazing() on phy per phy func type
    ZooGrazeGrowth = zoogrowth_zoo(GrazingMatrix, P, pfn ,zn)  # Zoo grazing * P summed up for zoo func types

    Losses = [p[i].losses(int_MLD, PhytoGrazed[i], DiffusiveMixing) for i in range(pfn)]

    # Zooplankton Growth & Grazing
    ZooMixing = [Z[i] * ActiveMixing for i in range(zn)] # Zooplankton actively stay within the MLD


    ZooGrowth = [z[i].zoogrowth(ZooGrazeGrowth[i]) for i in range(zn)]
    UnassimilatedProduction = [z[i].unassimilatedfeeding(ZooGrazeGrowth[i]) for i in range(zn)]


    y0 = y0 - sum([P[i ] *Gains[i] for i in range(pfn)]) # Nitrate drawdown

    y1 = y1 - sum([p[i].silicatedrawdown(P[i] ,Gains[i]) for i in range(pfn)]) # Silicate drawdown

    y2 = y2 + sum(UnassimilatedProduction) + sum([p[i].mortality(P[i]) for i in range(pfn)]) # Detritus


    zoo = [ZooGrowth[i ] -ZooMortality[i ] -ZooMixing[i] for i in range(zn)]    # Zooplankton losses due to mortality and mixing
    phy = [P[i] * (Gains[i] - Losses[i]) for i in range(pfn)]  # Phytoplankton growth

    out = [y0, y1, y2] + zoo + phy
    return out

