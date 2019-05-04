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
MLDfile = pandas.read_csv('Forcing/MLD2015_R1.csv')
NOXfile = pandas.read_csv('Forcing/NO2NO3_R1.csv')
SiOXfile = pandas.read_csv('Forcing/SiOH_R1.csv')
SSTfile = pandas.read_csv('Forcing/SST_R1.csv')
PARfile = pandas.read_csv('Forcing/PAR_R1.csv')

MLD_monthly_median = MLDfile.groupby('month').mean()
MLD = list(MLD_monthly_median['MLD'])
MLD.append(MLD[0])

NOX_monthly_median = NOXfile.groupby('month').mean()
NOX = list(NOX_monthly_median['NO2NO3'])
NOX.append(NOX[0])

SiOX_monthly_median = SiOXfile.groupby('month').mean()
SiOX = list(SiOX_monthly_median['SiOH'])
SiOX.append(SiOX[0])

SST_monthly_median = SSTfile.groupby('month').mean()
SST = list(SST_monthly_median['SST'])
SST.append(SST[0])

PAR_monthly_median = PARfile.groupby('month').mean()
PAR = list(PAR_monthly_median['PAR'])
PAR.append(PAR[0])

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
    if int_SIOX < 0: int_SIOX = 0  # do not allow negative Silicate values
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
    # SiRemineralization = paras['deltaD_Si'].value * D

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
    Rj = [z[i].ressourcedensity(P) for i in range(zn)] # total available ressource density for Zooplankton i

    Itot = [z[i].itot(Rj[i]) for i in range(zn)]

    PhytoGrazed = [p[i].zoograzing(Itot, Rj, P[i], Z) for i in range(pfn)] #returns phyto grazed per type

    ZooMixing = [Z[i] * K_Z for i in range(zn)]
    ZooMortality = [z[i].zoomortality(Z[i]) for i in range(zn)] #THIS IS IT!

    AssimilatedGrazing = [z[i].assimgrazing(Itot[i], Z[i]) for i in range(zn)]
    UnassimilatedGrazing = [z[i].unassimilatedgrazing(Itot[i], Z[i]) for i in range(zn)]

    InterZooPredation = [z[i].interzoograze(i, Z) for i in range(zn)]
    HigherOrderPredation = [z[i].higherorderpred(Z[i]) for i in range(zn)]

    # Phytoplankton losses
    PhytoMortality = [p[i].mortality(P[i]) for i in range(pfn)]
    PhytoSinking = [p[i].sinking(int_MLD, P[i]) for i in range(pfn)]
    PhytoMixing = [P[i] * K for i in range(pfn)]

    y0 = NRemineralization + NMixing - sum(Gains)  # Nitrate drawdown
    y1 = SiMixing - sum(SilicateDrawdown)  # Silicate drawdown
    y2 = sum(UnassimilatedGrazing) + sum(ZooMortality) + sum(PhytoMortality) - NRemineralization - DetritusMixing   # Detritus

    zoo = [AssimilatedGrazing[i] + InterZooPredation[i] - ZooMixing[i] - ZooMortality[i] - HigherOrderPredation[i] for i in range(zn)]   # Zooplankton losses due to mortality and mixing
    phy = [Gains[i] - PhytoGrazed[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]  # Phytoplankton growth

    outputlist[0] = NRemineralization       - outputlist[0]
    outputlist[1] = NMixing                 - outputlist[1]
    outputlist[2] = 0      - outputlist[2]
    outputlist[3] = SiMixing                - outputlist[3]
    outputlist[4] = DetritusMixing          - outputlist[4]

    outputlist[5] = sum(Gains)              - outputlist[5]
    outputlist[6] = sum(SilicateDrawdown)   - outputlist[6]
    outputlist[7] = sum(PhytoMortality)     - outputlist[7]
    outputlist[8] = sum(PhytoGrazed)        - outputlist[8]
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