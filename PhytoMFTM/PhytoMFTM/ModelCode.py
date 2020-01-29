#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

def phytomftm_extendedoutput_forcing(x, t, paras, pClass, zClass, forcing):
    N = x[0]  # Nitrate
    Si = x[1]  # Silicate
    D = x[2]  # Detritus

    zn = paras['zoo_num'].value
    Z = [x[3 + j] for j in range(zn)]  # Zooplankton

    z = zClass

    pfn = paras['pfun_num'].value
    P = [x[3 + zn + i] for i in range(pfn)]  # Phytoplankton

    p = pClass

    outputlist = [x[3+zn+pfn+i] for i in range(20)]

    # Interpolations of Forcings
    int_MLD = 100 #forcing.MLD.return_interpvalattime(t)
    int_NOX = forcing.NOX.return_interpvalattime(t)
    #if int_NOX < 0. : int_NOX = 0.  # do not allow negative Nitrate values
    int_SIOX = forcing.SiOX.return_interpvalattime(t)
    #if int_SIOX < 0.: int_SIOX = 0.  # do not allow negative Silicate values
    int_PAR = forcing.PAR.return_interpvalattime(t)
    int_SST = forcing.SST.return_interpvalattime(t)
    #print(int_MLD,int_NOX,int_SST)
    # Derivatives of Forcings
    deriv_MLD = forcing.MLD.return_derivattime(t)

    int_X21 = forcing.X21.return_interpvalattime(t)
    deriv_X21 = -forcing.X21.return_derivattime(t)

    # Non-Phytoplankton related processes
    # Mixing Processes
    if forcing.type == 'MLD':
        #K = (paras['kappa'].value + max(deriv_MLD, 0)) / int_MLD  # i.e. there is constant mixing & increased loss with MLD shallowing
        K = paras['kappa'].value / int_X21 * 10
        K_Z = 0 #deriv_MLD / int_MLD  # i.e. concentration varies with the MLD depth
        U = 0 #i.e. there is no mixing via U, the upwelling parameter for box models
        Sink = 0

    elif forcing.type == 'box':
        K = 0  # i.e. there is no modification of mixing due to K
        K_Z = 0  # since there is no change in modeled depth, no losses due to MLD changes
        U = (paras['kappa'].value + max(deriv_X21,0)) / 100
        # #(paras['kappa'].value + max(deriv_MLD, 0)) / int_MLD  #  max(deriv_MLD, 0)) / int_MLD # upwelling according to MLD
        Sink = 0.005  # 0.5 m per day


    elif forcing.type == 'box_constantKappa':
        K = 0  # i.e. there is no modification of mixing due to K
        K_Z = 0  # since there is no change in modeled depth, no losses due to MLD changes
        U = 0.02 # max(deriv_MLD, 0)) / int_MLD # upwelling according to MLD
        Sink = 0
    elif forcing.type == 'box_stochasticKappa':
        import random
        stochastic = random.uniform(0, 0.03)
        K = 0  # i.e. there is no modification of mixing due to K
        K_Z = 0  # since there is no change in modeled depth, no losses due to MLD changes
        U = stochastic  # max(deriv_MLD, 0)) / int_MLD # upwelling according to MLD
        Sink = 0
    elif forcing.type == 'box_MLD_stochastic':
        import random
        stochastic = random.uniform(-0.01, 0.01)
        K = 0  # i.e. there is no modification of mixing due to K
        K_Z = 0  # since there is no change in modeled depth, no losses due to MLD changes
        U = max((paras['kappa'].value + max(deriv_MLD,
                                        0)) / int_MLD + stochastic, 0) # max(deriv_MLD, 0)) / int_MLD # upwelling according to MLD
        Sink = 0
    else:
        raise('wrong forcing.type in forcing class, check forcing call')

    # Remineralisation
    NRemineralization = paras['deltaD_N'].value * D
    # SiRemineralization = paras['deltaD_Si'].value * D

    # Detritus
    DetritusMixing = D * max(K, Sink)

    # Nutrient Mixing
    NMixing = max(K, U) * (int_NOX - N)
    SiMixing = max(K, U) * (int_SIOX - Si)

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]

    # Light and Temperature
    LightHarvesting = [p[i].lightharvesting(int_MLD, int_PAR) for i in range(pfn)]
    # LightHarvesting = [p[i].smithpi(int_MLD, int_PAR, P) for i in range(pfn)]
    TemperatureDepGrowth = [p[i].tempdepgrowth(int_SST) for i in range(pfn)]
    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth[i], P[i]) for i in range(pfn)]
    SilicateDrawdown = [p[i].silicatedrawdown(Gains[i]) for i in range(pfn)]

    # Zooplankton Grazing:
    Gj = [z[j].zoofeeding(P, Z, func='anderson') for j in range(zn)]  # feeding probability for all food
    PhytoGrazed = [p[i].zoograzing(Gj, P[i], Z) for i in range(pfn)]  # returns phyto grazed per type

    ZooMixing = [Z[i] * K_Z for i in range(zn)]
    ZooMortality = [z[i].zoomortality(Z[i]) for i in range(zn)]  # THIS IS IT!

    ZooFeeding = [z[j].fullgrazing(Gj[j], P, Z, Z[j]) for j in range(zn)]

    AssimilatedGrazing = [z[j].assimgrazing(ZooFeeding[j]) for j in range(zn)]
    UnassimilatedGrazing = [z[j].unassimilatedgrazing(ZooFeeding[j]) for j in range(zn)]

    InterZooPredation = [z[j].interzoograze(Gj, Z, Z[j]) for j in range(zn)]
    HigherOrderPredation = [z[j].higherorderpred(Z[j]) for j in range(zn)]

    # Phytoplankton losses
    PhytoMortality = [p[i].mortality(P[i]) for i in range(pfn)]
    PhytoSinking = [p[i].sinking(int_MLD, P[i]) for i in range(pfn)]
    PhytoMixing = [P[i] * max(K, Sink) for i in range(pfn)]

    y0 = NRemineralization + NMixing - sum(Gains)  # Nitrate draw down
    y1 = SiMixing - sum(SilicateDrawdown)  # Silicate draw down
    y2 = sum(UnassimilatedGrazing) + sum(ZooMortality) + sum(PhytoMortality) - NRemineralization - DetritusMixing   # Detritus

    phy = [Gains[i] - PhytoGrazed[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]  # Phytoplankton growth
    zoo = [AssimilatedGrazing[j] - InterZooPredation[j] - ZooMixing[j] - ZooMortality[j] - HigherOrderPredation[j] for j in range(zn)]   # Zooplankton losses due to mortality and mixing

    #print(phy)
    outputlist[0] = U                 - outputlist[0]
    outputlist[1] = NMixing                 - outputlist[1]
    outputlist[2] = K                       - outputlist[2]
    outputlist[3] = sum(Gains)              - outputlist[3]
    outputlist[4] = U                       - outputlist[4]

    outputlist[5] = int_NOX                 - outputlist[5]
    outputlist[6] = N                       - outputlist[6]
    outputlist[7] = sum(P)                  - outputlist[7]
    outputlist[8] = sum(PhytoGrazed)        - outputlist[8]
    outputlist[9] = sum(PhytoMixing)        - outputlist[9]
    outputlist[10] = sum(PhytoSinking)      - outputlist[10]

    outputlist[11] = sum(Z)                 - outputlist[11]
    outputlist[12] = sum(AssimilatedGrazing)- outputlist[12]
    outputlist[13] = sum(InterZooPredation) - outputlist[13]
    outputlist[14] = sum(ZooMixing)         - outputlist[14]
    outputlist[15] = sum(ZooMortality)      - outputlist[15]
    outputlist[16] = sum(HigherOrderPredation)   - outputlist[16]

    outputlist[17] = sum(PhytoMortality)    - outputlist[17]
    outputlist[18] = deriv_X21 - outputlist[18]
    outputlist[19] = int_X21 - outputlist[19]

    out = [y0, y1, y2] + zoo + phy + outputlist
    return np.array(out)



def phytomftm_chemostat(x, t, paras, pClass, zClass, forcing):
    N = x[0]  # Nutrient 1
    Si = x[1]  # Nutrient 2

    D = x[2]  # Detritus

    zn = paras['zoo_num'].value
    Z = [x[3 + j] for j in range(zn)]  # Zooplankton

    z = zClass

    pfn = paras['pfun_num'].value
    P = [x[3 + zn + i] for i in range(pfn)]  # Phytoplankton

    p = pClass

    outputlist = [x[3+zn+pfn+i] for i in range(20)]

    int_MLD = 1
    deriv_MLD = 0
    # Interpolations of Forcings
    int_NOX = 1
    #if int_NOX < 0. : int_NOX = 0.  # do not allow negative Nitrate values
    int_SIOX = 0
    #if int_SIOX < 0.: int_SIOX = 0.  # do not allow negative Silicate values
    int_PAR = 40
    int_SST = 25

    # Non-Phytoplankton related processes
    # "Physics"
    if forcing.type == 'flowthrough':
        flow = paras['flow'].value
    elif forcing.type == 'batch':
        flow = 0
    else:
        raise('wrong forcing.type in forcing class, check forcing call')

    # Remineralisation
    NRemineralization = paras['deltaD_N'].value * D

    # Detritus
    DetritusMixing = D * flow

    # Nutrient Mixing
    NInput = int_NOX * flow
    SiInput = int_SIOX * flow

    # Phytoplankton related processes
    # Nutrient uptake
    N_Uptake = [p[i].n_uptake(N) for i in range(pfn)]
    Si_Uptake = [p[i].si_uptake(Si) for i in range(pfn)]

    # Light and Temperature
    # LightHarvesting = [p[i].lightharvesting(int_MLD, int_PAR) for i in range(pfn)]
    LightHarvesting = [p[i].smithpi(int_MLD, int_PAR, P) for i in range(pfn)]
    TemperatureDepGrowth = [p[i].tempdepgrowth(int_SST) for i in range(pfn)]
    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth[i], P[i]) for i in range(pfn)]
    SilicateDrawdown = [p[i].silicatedrawdown(Gains[i]) for i in range(pfn)]

    # Zooplankton Grazing:
    Gj = [z[j].zoofeeding(P, Z, func='anderson') for j in range(zn)]  # feeding probability for all food
    PhytoGrazed = [p[i].zoograzing(Gj, P[i], Z) for i in range(pfn)]  # returns phyto grazed per type

    ZooMixing = [Z[i] * flow for i in range(zn)]
    ZooMortality = [z[i].zoomortality(Z[i]) for i in range(zn)]  # THIS IS IT!

    ZooFeeding = [z[j].fullgrazing(Gj[j], P, Z, Z[j]) for j in range(zn)]

    AssimilatedGrazing = [z[j].assimgrazing(ZooFeeding[j]) for j in range(zn)]
    UnassimilatedGrazing = [z[j].unassimilatedgrazing(ZooFeeding[j]) for j in range(zn)]

    InterZooPredation = [z[j].interzoograze(Gj, Z, Z[j]) for j in range(zn)]
    HigherOrderPredation = [z[j].higherorderpred(Z[j]) for j in range(zn)]

    # Phytoplankton losses
    PhytoMortality = [p[i].mortality(P[i]) for i in range(pfn)]
    PhytoMixing = [P[i] * flow for i in range(pfn)]

    y0 = NRemineralization + NInput - sum(Gains)  # Nitrate draw down
    y1 = SiInput - sum(SilicateDrawdown)  # Silicate draw down
    y2 = sum(UnassimilatedGrazing) + sum(ZooMortality) + sum(PhytoMortality) - NRemineralization - DetritusMixing   # Detritus

    phy = [Gains[i] - PhytoGrazed[i] - PhytoMortality[i] - PhytoMixing[i] for i in range(pfn)]  # Phytoplankton growth
    zoo = [AssimilatedGrazing[j] - InterZooPredation[j] - ZooMixing[j] - ZooMortality[j] - HigherOrderPredation[j] for j in range(zn)]   # Zooplankton losses due to mortality and mixing

    #print(phy)
    outputlist[0] = int_MLD                 - outputlist[0]
    outputlist[1] = NInput                 - outputlist[1]
    outputlist[2] = flow                       - outputlist[2]
    outputlist[3] = sum(Gains)              - outputlist[3]
    outputlist[4] = 0                       - outputlist[4]

    outputlist[5] = int_NOX                 - outputlist[5]
    outputlist[6] = N                       - outputlist[6]
    outputlist[7] = sum(P)                  - outputlist[7]
    outputlist[8] = sum(PhytoGrazed)        - outputlist[8]
    outputlist[9] = sum(PhytoMixing)        - outputlist[9]
    outputlist[10] = sum(P)      - outputlist[10]

    outputlist[11] = sum(Z)                 - outputlist[11]
    outputlist[12] = sum(AssimilatedGrazing)- outputlist[12]
    outputlist[13] = sum(InterZooPredation) - outputlist[13]
    outputlist[14] = sum(ZooMixing)         - outputlist[14]
    outputlist[15] = sum(ZooMortality)      - outputlist[15]
    outputlist[16] = sum(HigherOrderPredation)   - outputlist[16]

    outputlist[17] = sum(PhytoMortality)    - outputlist[17]
    outputlist[18] = 0 - outputlist[18]
    outputlist[19] = 0 - outputlist[19]

    out = [y0, y1, y2] + zoo + phy + outputlist
    return np.array(out)



def phytomftm_extendedoutput(x, t, paras, pClass, zClass):
    N = x[0]  # Nitrate
    Si = x[1]  # Silicate
    D = x[2]  # Detritus

    zn = paras['zoo_num'].value
    Z = [x[3 + j] for j in range(zn)]  # Zooplankton

    z = zClass

    pfn = paras['pfun_num'].value
    P = [x[3 + zn + i] for i in range(pfn)]  # Phytoplankton

    p = pClass

    outputlist = [x[3+zn+pfn+i] for i in range(20)]

    # Interpolations of Forcings
    int_MLD = dailyinterp(MLD, t, kind=kindmld, k=kmld, s=smld)
    int_NOX = dailyinterp(NOX, t, kind=kindn0x, k=kn0x, s=sn0x)
    #if int_NOX < 0. : int_NOX = 0.  # do not allow negative Nitrate values
    int_SIOX = dailyinterp(SiOX, t, kind=kindsi0x, k=ksi0x, s=ssi0x)
    #if int_SIOX < 0.: int_SIOX = 0.  # do not allow negative Silicate values
    int_PAR = dailyinterp(PAR, t, kind=kindpar, k=kpar, s=spar)
    int_SST = dailyinterp(SST, t,  kind=kindsst, k=ksst, s=ssst)

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
    # LightHarvesting = [p[i].lightharvesting(int_MLD, int_PAR) for i in range(pfn)]
    LightHarvesting = [p[i].smithpi(int_MLD, int_PAR, P) for i in range(pfn)]
    TemperatureDepGrowth = [p[i].tempdepgrowth(int_SST) for i in range(pfn)]
    # Phytoplankton Growth
    Gains = [p[i].gains(N_Uptake[i], Si_Uptake[i], LightHarvesting[i], TemperatureDepGrowth[i], P[i]) for i in range(pfn)]
    SilicateDrawdown = [p[i].silicatedrawdown(Gains[i]) for i in range(pfn)]

    # Zooplankton Grazing:
    Gj = [z[j].zoofeeding(P, Z, func='anderson') for j in range(zn)]  # feeding probability for all food
    PhytoGrazed = [p[i].zoograzing(Gj, P[i], Z) for i in range(pfn)]  # returns phyto grazed per type

    ZooMixing = [Z[i] * K_Z for i in range(zn)]
    ZooMortality = [z[i].zoomortality(Z[i]) for i in range(zn)]  # THIS IS IT!

    ZooFeeding = [z[j].fullgrazing(Gj[j], P, Z, Z[j]) for j in range(zn)]

    AssimilatedGrazing = [z[j].assimgrazing(ZooFeeding[j]) for j in range(zn)]
    UnassimilatedGrazing = [z[j].unassimilatedgrazing(ZooFeeding[j]) for j in range(zn)]

    InterZooPredation = [z[j].interzoograze(Gj, Z, Z[j]) for j in range(zn)]
    HigherOrderPredation = [z[j].higherorderpred(Z[j]) for j in range(zn)]

    # Phytoplankton losses
    PhytoMortality = [p[i].mortality(P[i]) for i in range(pfn)]
    PhytoSinking = [p[i].sinking(int_MLD, P[i]) for i in range(pfn)]
    PhytoMixing = [P[i] * K for i in range(pfn)]

    y0 = NRemineralization + NMixing - sum(Gains)  # Nitrate draw down
    y1 = SiMixing - sum(SilicateDrawdown)  # Silicate draw down
    y2 = sum(UnassimilatedGrazing) + sum(ZooMortality) + sum(PhytoMortality) - NRemineralization - DetritusMixing   # Detritus

    phy = [Gains[i] - PhytoGrazed[i] - PhytoMortality[i] - PhytoMixing[i] - PhytoSinking[i] for i in range(pfn)]  # Phytoplankton growth
    zoo = [AssimilatedGrazing[j] - InterZooPredation[j] - ZooMixing[j] - ZooMortality[j] - HigherOrderPredation[j] for j in range(zn)]   # Zooplankton losses due to mortality and mixing

    #print(phy)
    outputlist[0] = NRemineralization       - outputlist[0]
    outputlist[1] = NMixing                 - outputlist[1]
    outputlist[2] = 1                       - outputlist[2] #K
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
    outputlist[17] = sum(InterZooPredation) - outputlist[17]
    outputlist[18] = 0 - outputlist[18]
    outputlist[19] = sum(HigherOrderPredation) - outputlist[19]

    out = [y0, y1, y2] + zoo + phy + outputlist
    return np.array(out)