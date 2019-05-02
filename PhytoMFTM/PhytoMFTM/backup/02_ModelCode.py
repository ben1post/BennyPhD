#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import pandas
from PhytoMFTM.AuxFuncs import dailyinterp, firstderivspl



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



def simpleN2P2ZD_extendedoutput(x, t, paras, pClass, zClass):
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
    GrazingMatrix = grazing(Ptot, Z, pfn,zn, z) # return all zooplankton feeding on all phytoplankton types

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

    zoo = [ZooGrowth[i ] -ZooMortality[i ] -ZooMixing[i] for i in range(zn)]

    phy = [P[i] * (Gains[i] - Losses[i]) for i in range(pfn)]  # Phytoplankton growth

    outputlist = [DiffusiveMixing, ActiveMixing, NRemineralization, SiRemineralization,
                  DetritusMixing, NMixing, SiMixing, sum(ZooMortality), sum(N_Uptake), sum(Si_Uptake),
                  sum(LightHarvesting), TemperatureDepGrowth, sum(Gains), sum(PhytoGrazed), sum(ZooGrazeGrowth),
                  sum(Losses), sum(ZooMixing), sum(ZooGrowth), sum(UnassimilatedProduction)]

    out = [y0, y1, y2] + zoo + phy + outputlist
    return out



def simpleN2P2ZD_extendedoutput_constantinput(x, t, paras, pClass, zClass):
    N = x[0]  # Nitrate
    Si = x[1]  # Silicate
    D = x[2]  # Detritus

    zn = paras['zoo_num'].value
    Z = [x[3 +i] for i in range(zn)]  # Zooplankton

    z = zClass

    pfn = paras['pfun_num'].value
    P = [x[3 + zn + i] for i in range(pfn)]  # Phytoplankton

    p = pClass

    outputlist = [x[3+zn+pfn+i] for i in range(17)]

    # Interpolations of Forcings
    int_MLD = dailyinterp(MLD, t, kind=kindmld, k=kmld, s=smld)
    int_NOX = dailyinterp(NOX, t, kind=kindn0x, k=kn0x, s=sn0x)
    if int_NOX < 0 : int_NOX = 0  # do not allow negative Nitrate values
    int_SIOX = dailyinterp(SiOX, t, kind=kindsi0x, k=ksi0x, s=ssi0x)
    int_PAR = dailyinterp(PAR, t, kind=kindpar, k=kpar, s=spar)
    int_SST = dailyinterp(SST, t, kind=kindsst, k=ksst, s=ssst)

    # Derivatives of Forcings
    deriv_MLD = firstderivspl(MLD, t, k=kmld, s=smld)

    # Non-Phytoplankton related processes
    # Mixing Processes
    DiffusiveMixing = (paras['kappa'].value + max(deriv_MLD, 0)) / int_MLD  # i.e. there is constant mixing & increased loss with MLD shallowing
    ActiveMixing = deriv_MLD / int_MLD  # i.e. concentration varies with the MLD depth

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

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]
    # Light and Temperature
    LightHarvesting = [p[i].lightharvesting(int_MLD, int_PAR) for i in range(pfn)]
    TemperatureDepGrowth = [p[i].tempdepgrowth(int_SST) for i in range(pfn)]
    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth[i], P[i]) for i in range(pfn)]
    SilicateDrawdown = [p[i].silicatedrawdown(Gains[i]) for i in range(pfn)]

    # Zooplankton Grazing:
    GrazeZoo = [z[i].zoograzing(P) for i in range(zn)]
    ZooItots = [x[0] for x in GrazeZoo]
    RperZ = [x[1] for x in GrazeZoo]
    # Zooplankton Growth & Grazing
    ZooMixing = [Z[i] * ActiveMixing for i in range(zn)]
    ZooGrowth = [z[i].zoogrowth(ZooItots[i], Z[i]) for i in range(zn)]
    UnassimilatedProduction = [z[i].unassimilatedfeeding(ZooItots[i]) for i in range(zn)]

    # Phytoplankon being grazed
    PhytoGrazed = [p[i].grazedphyto(Itot=ZooItots, P=P[i], R=RperZ) for i in range(pfn)]  # return intake rate per type
    # Phytoplankton losses
    PhytoMortality = [p[i].mortality(P[i]) for i in range(pfn)]
    PhytoSinking = [p[i].sinking(int_MLD, P[i]) for i in range(pfn)]
    PhytoMixing = [p[i].mixing(DiffusiveMixing, P[i]) for i in range(pfn)]
    PhytoGrazing = [p[i].grazed(PhytoGrazed[i], P[i]) for i in range(pfn)]  # returns biomass grazed by zoo per type

    y0 = NRemineralization + NMixing - sum(Gains)  # Nitrate drawdown
    y1 = SiRemineralization + SiMixing - sum(SilicateDrawdown)  # Silicate drawdown
    y2 = sum(ZooMortality) + sum(UnassimilatedProduction) + sum(PhytoMortality) - NRemineralization - SiRemineralization - DetritusMixing   # Detritus

    zoo = [ZooGrowth[i] - ZooMortality[i] - ZooMixing[i] for i in range(zn)]    # Zooplankton losses due to mortality and mixing
    phy = [Gains[i] - PhytoGrazing[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]  # Phytoplankton growth

    outputlist[0] = NRemineralization       - outputlist[0]
    outputlist[1] = NMixing                 - outputlist[1]
    outputlist[2] = SiRemineralization      - outputlist[2]
    outputlist[3] = SiMixing                - outputlist[3]
    outputlist[4] = DetritusMixing          - outputlist[4]

    outputlist[5] = sum(Gains)              - outputlist[5]
    outputlist[6] = sum(SilicateDrawdown)   - outputlist[6]
    outputlist[7] = sum(PhytoMortality)     - outputlist[7]
    outputlist[8] = sum(PhytoGrazing)       - outputlist[8]
    outputlist[9] = sum(PhytoMixing)        - outputlist[9]
    outputlist[10] = sum(PhytoSinking)      - outputlist[10]

    outputlist[11] = sum([LightHarvesting[i] * P[i] for i in range(pfn)])           - outputlist[11]
    outputlist[12] = sum([TemperatureDepGrowth[i] * P[i] for i in range(pfn)])      - outputlist[12]

    outputlist[13] = sum(ZooGrowth)         - outputlist[13]
    outputlist[14] = sum(ZooMortality)      - outputlist[14]
    outputlist[15] = sum(ZooMixing)         - outputlist[15]
    outputlist[16] = sum(UnassimilatedProduction)   - outputlist[16]

    out = [y0, y1, y2] + zoo + phy + outputlist
    return np.array(out)

def phytomftm_for(x, t, paras, pClass, zClass):
    N = x[0]  # Nitrate
    Si = x[1]  # Silicate
    D = x[2]  # Detritus

    z = zClass
    p = pClass

    zn = paras['zoo_num'].value
    pfn = paras['pfun_num'].value

    Z = [x[3 +i] for i in range(zn)]  # Zooplankton
    P = [x[3 + zn + i] for i in range(pfn)]  # Phytoplankton

    phy = [0 for i in range(pfn)]
    zoo = [0 for i in range(zn)]

    #outputlist = [x[3+zn+pfn+i] for i in range(17)]

    # Interpolations of Forcings
    int_MLD = dailyinterp(MLD, t, kind=kindmld, k=kmld, s=smld)
    int_NOX = dailyinterp(NOX, t, kind=kindn0x, k=kn0x, s=sn0x)
    if int_NOX < 0 : int_NOX = 0  # do not allow negative Nitrate values
    int_SIOX = dailyinterp(SiOX, t, kind=kindsi0x, k=ksi0x, s=ssi0x)
    int_PAR = dailyinterp(PAR, t, kind=kindpar, k=kpar, s=spar)
    int_SST = dailyinterp(SST, t, kind=kindsst, k=ksst, s=ssst)

    # Derivatives of Forcings
    deriv_MLD = firstderivspl(MLD, t, k=kmld, s=smld)

    # Non-Phytoplankton related processes
    # Mixing Processes
    DiffusiveMixing = (paras['kappa'].value + max(deriv_MLD, 0)) / int_MLD  # i.e. there is constant mixing & increased loss with MLD shallowing
    ActiveMixing = deriv_MLD / int_MLD  # i.e. concentration varies with the MLD depth

    # Remineralisation
    NRemineralization = paras['deltaD_N'].value * D
    SiRemineralization = paras['deltaD_Si'].value * D

    # Detritus
    DetritusMixing = D * DiffusiveMixing

    # Nutrient Mixing
    NMixing = DiffusiveMixing * (int_NOX - N)
    SiMixing = DiffusiveMixing * (int_SIOX - Si)

    y0 = NRemineralization + NMixing
    y1 = SiRemineralization + SiMixing
    y2 = - NRemineralization - SiRemineralization - DetritusMixing  # Detritus

    ####################################################################################################################

    for i in range(pfn):
        # Nutrient uptake
        N_Uptake = p[i].n_uptake(N)
        Si_Uptake = p[i].si_uptake(Si)

        # Light and Temperature
        LightHarvesting = p[i].lightharvesting(int_MLD, int_PAR)
        TemperatureDepGrowth = p[i].tempdepgrowth(int_SST)

        # Phytoplankton Growth
        Gains = p[i].gains(N_Uptake, Si_Uptake, LightHarvesting, TemperatureDepGrowth, P[i])
        SilicateDrawdown = p[i].silicatedrawdown(Gains)

        # Phytoplankton losses
        PhytoMortality = p[i].mortality(P[i])
        PhytoSinking = p[i].sinking(int_MLD, P[i])
        PhytoMixing = p[i].mixing(DiffusiveMixing, P[i])

        phy[i] = (Gains - PhytoMortality - PhytoMixing - PhytoSinking)

        y0 = y0 - Gains  # Nitrate drawdown
        y1 = y1 - SilicateDrawdown  # Silicate drawdown
        y2 = y2 + PhytoMortality

    for i in range(zn):
        # Zooplankton
        ZooMortality = z[i].zoomortality(Z[i])
        # Zooplankton Growth & Grazing
        ZooMixing = Z[i] * ActiveMixing

        # Zooplankton Grazing:
        GrazeZoo = z[i].zoograzing(P)
        ZooItots = GrazeZoo[0]
        RperZ = GrazeZoo[1]
        # Zooplankton Growth & Grazing
        ZooGrowth = z[i].zoogrowth(ZooItots, Z[i])
        UnassimilatedProduction = z[i].unassimilatedfeeding(ZooItots)

        for j in range(pfn):
            # Phytoplankon being grazed
            PhytoGrazed = p[j].grazedphyto_for(Itot=ZooItots, P=P[j], R=RperZ)
            PhytoGrazing = p[j].grazed(PhytoGrazed, P[j])
            phy[j] = phy[j] - PhytoGrazing

        y2 = y2 + ZooMortality + UnassimilatedProduction
        zoo[i] = ZooGrowth - ZooMortality - ZooMixing

    out = [y0, y1, y2] + zoo + phy
    return np.array(out)


def OLD_phytomftm_extendedoutput(x, t, paras, pClass, zClass):
    N = x[0]  # Nitrate
    Si = x[1]  # Silicate
    D = x[2]  # Detritus

    zn = paras['zoo_num'].value
    Z = [0 for i in range(zn)]  # [x[3 +i] for i in range(zn)]  # Zooplankton

    z = zClass

    pfn = paras['pfun_num'].value
    P = [x[3 + zn + i] for i in range(pfn)]  # Phytoplankton

    p = pClass

    outputlist = [x[3+zn+pfn+i] for i in range(17)]

    # Interpolations of Forcings
    int_MLD = dailyinterp(MLD, t, kind=kindmld, k=kmld, s=smld)
    int_NOX = dailyinterp(NOX, t, kind=kindn0x, k=kn0x, s=sn0x)
    if int_NOX < 0 : int_NOX = 0  # do not allow negative Nitrate values
    int_SIOX = dailyinterp(SiOX, t, kind=kindsi0x, k=ksi0x, s=ssi0x)
    int_PAR = dailyinterp(PAR, t, kind=kindpar, k=kpar, s=spar)
    int_SST = dailyinterp(SST, t, kind=kindsst, k=ksst, s=ssst)

    # Derivatives of Forcings
    deriv_MLD = firstderivspl(MLD, t, k=kmld, s=smld)

    # Non-Phytoplankton related processes
    # Mixing Processes
    DiffusiveMixing = (paras['kappa'].value + max(deriv_MLD, 0)) / int_MLD
    ZooMixing = deriv_MLD / int_MLD  # i.e. concentration varies with the MLD depth

    # Remineralisation
    NRemineralization = paras['deltaD_N'].value * D
    SiRemineralization = 0  # paras['deltaD_Si'].value * D

    # Detritus
    DetritusMixing = D * DiffusiveMixing

    # Nutrient Mixing
    NMixing = DiffusiveMixing * (int_NOX - N)
    SiMixing = DiffusiveMixing * (int_SIOX - Si)

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]




    # Light and Temperature
    LightHarvesting = [p[i].lightharvesting(int_MLD, int_PAR) for i in range(pfn)]
    TemperatureDepGrowth = [p[i].tempdepgrowth(int_SST) for i in range(pfn)]

    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth[i]) for i in range(pfn)]
    SilicateDrawdown = [p[i].silicatedrawdown(Gains[i]) for i in range(pfn)]

    # Zooplankton Grazing:
    RperZ = [z[i].ressourcedensity(P) for i in range(zn)]
    ZooGrazing = [z[i].zoograzing(RperZ[i], P) for i in range(zn)]

    # Zooplankton Growth & Grazing
    ZooMixing = [0 for i in range(zn)]  # [Z[i] * ActiveMixing for i in range(zn)]
    ZooGrowth = [z[i].zoogrowth(ZooGrazing[i]) for i in range(zn)]
    # Zooplankton
    ZooMortality = [z[i].zoomortality(Z[i]) for i in range(zn)]

    UnassimilatedProduction = [z[i].unassimilatedfeeding(ZooGrazing[i]) for i in range(zn)]

    # Phytoplankon being grazed
    PhytoGrazed = [0 for i in range(pfn)]  # [p[i].grazedphyto(ZooGrazing, P[i], RperZ) for i in range(pfn)]  # return intake rate per type
    # Phytoplankton losses
    PhytoMortality = [p[i].mortality() for i in range(pfn)]
    PhytoSinking = [p[i].sinking(int_MLD) for i in range(pfn)]
    PhytoMixing = [p[i].mixing(DiffusiveMixing) for i in range(pfn)]
    PhytoGrazing = [p[i].grazed(PhytoGrazed[i]) for i in range(pfn)]  # returns biomass grazed by zoo per type
    Losses = [-PhytoGrazing[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]

    y0 = NRemineralization + NMixing - sum([Gains[i] * P[i] for i in range(pfn)])  # Nitrate drawdown
    y1 = SiRemineralization + SiMixing - sum([SilicateDrawdown[i] * P[i] for i in range(pfn)])  # Silicate drawdown

    y2 = sum([PhytoMortality[i] * P[i] for i in range(pfn)]) - \
        NRemineralization - SiRemineralization - DetritusMixing + \
        sum([(UnassimilatedProduction[i] * Z[i]) + ZooMortality[i] for i in range(zn)])  # Detritus

    zoo = [0 for i in range(zn)]  # [(Z[i] * (ZooGrowth[i] - ZooMixing[i])) - ZooMortality[i] for i in range(zn)]    # Zooplankton losses due to mortality and mixing
    phy = [P[i] * (Gains[i] - Losses[i]) for i in range(pfn)]  # Phytoplankton growth

    outputlist[0] = NRemineralization       - outputlist[0]
    outputlist[1] = NMixing                 - outputlist[1]
    outputlist[2] = SiRemineralization      - outputlist[2]
    outputlist[3] = SiMixing                - outputlist[3]
    outputlist[4] = DetritusMixing          - outputlist[4]

    outputlist[5] = sum(Gains)              - outputlist[5]
    outputlist[6] = sum(SilicateDrawdown)   - outputlist[6]
    outputlist[7] = sum(PhytoMortality)     - outputlist[7]
    outputlist[8] = sum(PhytoGrazing)       - outputlist[8]
    outputlist[9] = sum(PhytoMixing)        - outputlist[9]
    outputlist[10] = sum(PhytoSinking)      - outputlist[10]

    outputlist[11] = sum([LightHarvesting[i] * P[i] for i in range(pfn)])           - outputlist[11]
    outputlist[12] = sum([TemperatureDepGrowth[i] * P[i] for i in range(pfn)])      - outputlist[12]

    outputlist[13] = sum(ZooGrowth)         - outputlist[13]
    outputlist[14] = sum(ZooMortality)      - outputlist[14]
    outputlist[15] = sum(ZooMixing)         - outputlist[15]
    outputlist[16] = sum(UnassimilatedProduction)   - outputlist[16]

    out = [y0, y1, y2] + zoo + phy + outputlist
    return np.array(out)


def phytomftm_extendedoutput(x, t, paras, pClass, zClass):
    N = x[0]  # Nitrate
    Si = x[1]  # Silicate
    D = x[2]  # Detritus

    zn = paras['zoo_num'].value
    Z = [x[3 +i] for i in range(zn)]  # Zooplankton

    z = zClass

    pfn = paras['pfun_num'].value
    P = [x[3 + zn + i] for i in range(pfn)]  # Phytoplankton

    p = pClass

    outputlist = [x[3+zn+pfn+i] for i in range(17)]

    # Interpolations of Forcings
    int_MLD = dailyinterp(MLD, t, kind=kindmld, k=kmld, s=smld)
    int_NOX = dailyinterp(NOX, t, kind=kindn0x, k=kn0x, s=sn0x)
    if int_NOX < 0 : int_NOX = 0  # do not allow negative Nitrate values
    int_SIOX = dailyinterp(SiOX, t, kind=kindsi0x, k=ksi0x, s=ssi0x)
    int_PAR = dailyinterp(PAR, t, kind=kindpar, k=kpar, s=spar)
    int_SST = dailyinterp(SST, t, kind=kindsst, k=ksst, s=ssst)

    # Derivatives of Forcings
    deriv_MLD = firstderivspl(MLD, t, k=kmld, s=smld)

    # Non-Phytoplankton related processes
    # Mixing Processes
    K = (paras['kappa'].value + max(deriv_MLD, 0)) / int_MLD  # i.e. there is constant mixing & increased loss with MLD shallowing
    K_Z = deriv_MLD / int_MLD  # i.e. concentration varies with the MLD depth

    # Remineralisation
    NRemineralization = paras['deltaD_N'].value * D
    SiRemineralization = paras['deltaD_Si'].value * D

    # Detritus
    DetritusMixing = D * K

    # Nutrient Mixing
    NMixing = K * (int_NOX - N)
    SiMixing = K * (int_SIOX - Si)

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]

    # Light and Temperature
    LightHarvesting = [p[i].lightharvesting(int_MLD, int_PAR) for i in range(pfn)]
    TemperatureDepGrowth = [p[i].tempdepgrowth(int_SST) for i in range(pfn)]
    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth[i], P[i]) for i in range(pfn)]
    SilicateDrawdown = [p[i].silicatedrawdown(Gains[i]) for i in range(pfn)]

    # Zooplankton Grazing:
    Rj = [z[i].ressourcedensity(P) for i in range(zn)] # total available ressource density for Zooplankton 1
    #Rj = [sum(Rji) for i in range(zn)]
    Itot = [z[i].itot(Rj[i]) for i in range(zn)]

    PhytoGrazed = [p[i].zoograzing(Itot, Rj, P[i], Z) for i in range(pfn)] #returns phyto grazed per type
    # p1 [grazed by z1, z2], p2 [grazed by z1,z2]
    ZooMixing = [Z[i] * K_Z for i in range(zn)]
    ZooMortality = [z[i].zoomortality(Z) for i in range(zn)] #THIS IS IT!

    AssimilatedGrazing = [z[i].assimgrazing(Itot[i], Z[i]) for i in range(zn)]
    UnassimilatedGrazing = [z[i].unassimilatedgrazing(Itot[i], Z[i]) for i in range(zn)]

    # Phytoplankton losses
    PhytoMortality = [p[i].mortality(P[i]) for i in range(pfn)]
    PhytoSinking = [p[i].sinking(int_MLD, P[i]) for i in range(pfn)]
    PhytoMixing = [P[i] * K for i in range(pfn)]

    y0 = NRemineralization + NMixing - sum(Gains)  # Nitrate drawdown
    y1 = SiRemineralization + SiMixing - sum(SilicateDrawdown)  # Silicate drawdown
    y2 = sum(UnassimilatedGrazing) + sum(ZooMortality) + sum(PhytoMortality) - NRemineralization - SiRemineralization - DetritusMixing   # Detritus

    zoo = [AssimilatedGrazing[i] - ZooMortality[i] - ZooMixing[i] for i in range(zn)]    # Zooplankton losses due to mortality and mixing
    phy = [Gains[i] - PhytoGrazed[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]  # Phytoplankton growth

    outputlist[0] = NRemineralization       - outputlist[0]
    outputlist[1] = NMixing                 - outputlist[1]
    outputlist[2] = SiRemineralization      - outputlist[2]
    outputlist[3] = SiMixing                - outputlist[3]
    outputlist[4] = DetritusMixing          - outputlist[4]

    outputlist[5] = sum(Gains)              - outputlist[5]
    outputlist[6] = sum(SilicateDrawdown)   - outputlist[6]
    outputlist[7] = sum(PhytoMortality)     - outputlist[7]
    outputlist[8] = sum(PhytoGrazed)       - outputlist[8]
    outputlist[9] = sum(PhytoMixing)        - outputlist[9]
    outputlist[10] = sum(PhytoSinking)      - outputlist[10]

    outputlist[11] = sum([LightHarvesting[i] * P[i] for i in range(pfn)])           - outputlist[11]
    outputlist[12] = sum([TemperatureDepGrowth[i] * P[i] for i in range(pfn)])      - outputlist[12]

    outputlist[13] = sum(AssimilatedGrazing)         - outputlist[13]
    outputlist[14] = sum(ZooMortality)      - outputlist[14]
    outputlist[15] = sum(ZooMixing)         - outputlist[15]
    outputlist[16] = sum(UnassimilatedGrazing)   - outputlist[16]

    out = [y0, y1, y2] + zoo + phy + outputlist
    return np.array(out)
