#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from runs.modelrun_aggTS_RegimeComp import timedays, outarray1, outarray2, ms1, ms2

# TODO:
#  - import necessary functions & objects

import matplotlib
params = {'font.size': 7}
matplotlib.rcParams.update(params)


numcols = 4
f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, numcols, sharex='col')  # , sharey='row')


plt.setp((ax1, ax2, ax3, ax4), xticks=[1,60,120,180,240,300,365])

from matplotlib.ticker import MaxNLocator
for axe in (ax1, ax2, ax3, ax4):
    for i in range(numcols):
        axe[i].get_yaxis().set_major_locator(MaxNLocator(nbins=4))
        axe[i].tick_params(top=True, right=True)


def plotstuff(ms,outarray, zuplot, muplot, regime):
    # PLOTTING
    timedays_ly = timedays[1:366]
    outarray_ly = outarray[1460:1825]

    # color vectors
    #colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', 'grey']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [2, 2.5, 4, 5.5]

    ax1[zuplot].set_title(regime + 'model out')

    dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    dpm_cumsum = np.cumsum(dayspermonth) - np.array(dayspermonth)/2 #- 15
    #print(timedays_ly)


    # Figure 1
    N = ms.physics.forcing.verif.NO2NO3
    Nmean = ms.physics.forcing.verif.returnMeanVerifPerMonth(N,'NO3_NO2_USF_Box')
    #print(N)
    # N
    N_Max = 20  # np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) + np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) * 0.1
    ax1[zuplot].scatter(N['yday'], N['NO3_NO2_USF_Box'], label='data', alpha=0.5)
    ax1[zuplot].plot(dpm_cumsum, Nmean, label='mean data', alpha=1)
    ax1[zuplot].plot(timedays_ly, outarray_ly[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
    ax1[0].set_ylabel('Nutrients \n' '[µM N]', multialignment='center')
    ax1[zuplot].set_ylim(0, N_Max)
    ax1[zuplot].legend(fontsize='x-small')


    ChlConv = True
    # Phyto
    CtoChla = 50  # g/g
    MolarMassC = 12.0107
    CtoNratioPhyto = 6.625
    muMolartoChlaconvfactor = CtoChla / MolarMassC / CtoNratioPhyto  # Chla as mg/m-3 to


    chla = ms.physics.forcing.verif.FluorChla

    chla2 = ms.physics.forcing.verif.HPLC

    ax2[zuplot].scatter(chla['yday'], chla['Chlorophyll_Box'] * muMolartoChlaconvfactor, label='Fluor Chla', alpha=0.5)
    ax2[zuplot].scatter(chla2['yday'], chla2['Tchla'] * muMolartoChlaconvfactor / 100, label='HPLC Chla', alpha=0.5)


    chlamean = ms.physics.forcing.verif.returnMeanVerifPerMonth(chla,'Chlorophyll_Box')

    chla2mean = ms.physics.forcing.verif.returnMeanVerifPerMonth(chla2,'Tchla')

    ax2[zuplot].plot(dpm_cumsum, np.array(chlamean) * muMolartoChlaconvfactor, label='mean FlChla', alpha=1)
    ax2[zuplot].plot(dpm_cumsum, np.array(chla2mean) * muMolartoChlaconvfactor / 100, label='mean HPLC', alpha=1)

    Pall = outarray_ly[:,1]
    P_Max = 2 # np.max(Pall) + 0.9 * np.max(Pall)

    ax2[zuplot].plot(timedays_ly, Pall, c=colors[4], lw=lws[1], label='Model')
    ax2[zuplot].legend(fontsize='x-small')
    ax2[0].set_ylabel('Phytoplankton \n' '[µM N]', multialignment='center')
    ax2[zuplot].set_ylim(0, P_Max)


    # mg dry weight per cubic meter to µM of N
    mggramstograms = 1/1000
    Cperdryweight = 0.32
    # Wiebe et al. 1975 : Carbon was 31-33% ofdryweight
    molarmassCarbon = 12.01 #grams per mole
    CtonNratioZoo = 5.625
    mgDWtomuMolarZOO = mggramstograms / Cperdryweight / molarmassCarbon / CtonNratioZoo * 1000 # µM

    ZBM = ms.physics.forcing.verif.padZoo
    #print(ZBM)
    ZBMmean = ms.physics.forcing.verif.returnMeanVerifPerMonth(ZBM, 'BIOMASS_200')

    #print(ZBMmean)
    ax3[zuplot].scatter(ZBM['yday'], ZBM['BIOMASS_200'] * muMolartoChlaconvfactor * mgDWtomuMolarZOO, label='200 µM', alpha=0.5)
    ax3[zuplot].plot(dpm_cumsum, np.array(ZBMmean) * muMolartoChlaconvfactor * mgDWtomuMolarZOO, label='mean data',
                        alpha=1)

    # Z
    Zall = outarray_ly[:,2]
    Z_Max = 1. #np.max(Zall) + 0.1 * np.max(Zall)

    ax3[zuplot].plot(timedays_ly, Zall, c=colors[4], lw=lws[1], label='model')
    ax3[0].set_ylabel('Zooplankton \n' '[µM N]', multialignment='center')
    ax3[zuplot].tick_params('y')
    ax3[zuplot].set_ylim(0, Z_Max)
    #ax4[i_plot].set_title('Zooplankton')
    ax3[zuplot].legend(fontsize='x-small')


    # convert PN in µg/L to µM of Detritus!
    molarmassNitrogen = 14.0067
    mugperlitertomuMolarPN = 1 / molarmassNitrogen  # g/L -> mol/L -> µM

    # TODO: add other verification variables! ensure proper converison..
    PN = ms.physics.forcing.verif.PN

    PNmean = ms.physics.forcing.verif.returnMeanVerifPerMonth(PN, 'PON_ug_kg_Box')

    #print(np.array(PNmean),dpm_cumsum)

    ax4[zuplot].plot(dpm_cumsum, np.array(PNmean) /14.0067 , label='mean data', alpha=0.5)
    ax4[zuplot].scatter(PN['yday'], PN['PON_ug_kg_Box'] /14.0067 , label='PON data', alpha=0.5)

    D_Max = 1.5 #np.max(outarray_ly[:, 3]) + 0.2 * np.max(outarray_ly[:, 3])
    # D
    ax4[zuplot].plot(timedays_ly, outarray_ly[:, 3], c=colors[1], lw=lws[0], alpha=alphas[0], label='model')
    ax4[0].set_ylabel('Detritus \n' '[µM N]', multialignment='center')
    ax4[zuplot].set_ylim(0, D_Max)
    ax4[zuplot].set_xlabel('Day in year')
    # Legend
    ax4[zuplot].legend(fontsize='x-small')


    ## PHYSICS ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    ax1[muplot].set_title(regime + 'forcing')

    ax4[muplot].set_xlabel('Day in year')

    NOX = ms.physics.forcing.NOX.return_interpvalattime(timedays_ly)
    NOXdat = ms.physics.forcing.NOX.forcingfile

    NN = ms.physics.forcing.NOX.rawforcing

    ax1[muplot].scatter(NN['yday'], NN['NO3_NO2_USF_AtDepth'] , label='data', alpha=0.5)

    #print(NOX)
    #print(NOXdat)
    ax1[muplot].plot(timedays_ly, NOX, c=colors[5], lw=lws[0], alpha=alphas[0])
    ax1[muplot].scatter(dpm_cumsum, NOXdat[0:12], c=colors[5])
    ax1[2].set_ylabel('$N_0$ \n' '[µM]', multialignment='center')
    ax1[muplot].set_ylim(0., N_Max)
    #ax1[muplot].invert_yaxis()


    MLDRAW = ms.physics.forcing.X258.rawforcing

    ax2[muplot].scatter(MLDRAW['yday'], MLDRAW['x258depth'], label='data', alpha=0.5)

    MLD = ms.physics.forcing.X258.return_interpvalattime(timedays_ly)
    MLDdat = ms.physics.forcing.X258.forcingfile
    MLD_max = 150 #np.max(MLDdat) + 0.1 * np.max(MLDdat)
    ax2[muplot].plot(timedays_ly, MLD, c=colors[5], lw=lws[0], alpha=alphas[0])
    ax2[muplot].scatter(dpm_cumsum, MLDdat[0:12], c=colors[5])
    ax2[2].set_ylabel('$25.8_{iso}$ \n' '[m]', multialignment='center')
    ax2[muplot].set_ylim(0, MLD_max) # 400 for biotrans, 100 for Papa
    ax2[muplot].invert_yaxis()


    PARAW = ms.physics.forcing.PAR.rawforcing
    #print(PARAW)
    ax3[muplot].scatter(PARAW['yday'], PARAW['PAR'], label='data', alpha=0.5)

    PAR = ms.physics.forcing.PAR.return_interpvalattime(timedays_ly)
    PARdat = ms.physics.forcing.PAR.forcingfile
    PAR_max = 70 #np.max(PAR) + 0.1 * np.max(PAR)
    ax3[muplot].plot(timedays_ly, PAR, c=colors[5], lw=lws[0], alpha=alphas[0])
    ax3[muplot].scatter(dpm_cumsum, PARdat[0:12], c=colors[5])
    ax3[2].set_ylabel('PAR \n' '[E $m^{−2}$ $s^{−1}$]', multialignment='center')
    ax3[muplot].set_ylim(0, PAR_max)
    # ax1[muplot].invert_yaxis()

    TmldRAW = ms.physics.forcing.SST.rawforcing
    #print(TmldRAW)
    ax4[muplot].scatter(TmldRAW['yday'], TmldRAW['Temperature_Box'], label='data', alpha=0.5)

    Tmld = ms.physics.forcing.SST.return_interpvalattime(timedays_ly)
    Tmlddat = ms.physics.forcing.SST.forcingfile
    Tmld_max = 30 #np.max(Tmld) + 0.1 * np.max(Tmld)
    ax4[muplot].plot(timedays_ly, Tmld, c=colors[5], lw=lws[0], alpha=alphas[0])
    ax4[muplot].scatter(dpm_cumsum, Tmlddat[0:12], c=colors[5])
    ax4[2].set_ylabel('$T_{MLD}$ \n' '[°C]', multialignment='center')
    ax4[muplot].set_ylim(18, Tmld_max)
    # ax1[muplot].invert_yaxis()



plotstuff(ms1, outarray1, 0, 2, 'Regime 1')

plotstuff(ms2, outarray2, 1, 3, 'Regime 2')

# Defining custom 'xlim' and 'ylim' values.
xlim = (0, 365)

# Setting the values for all axes.
plt.setp((ax1, ax2, ax3, ax4), xlim=xlim)

f1.align_ylabels()


#plt.subplots_adjust(vspace=0, hspace=0)

#plt.tight_layout(rect=[1, 1, 1, 1])
plt.show()