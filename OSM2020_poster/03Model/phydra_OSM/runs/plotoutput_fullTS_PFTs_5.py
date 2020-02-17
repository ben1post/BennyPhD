#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

import matplotlib
#matplotlib.use('WxAgg')

import matplotlib.pyplot as plt

from runs.modelrun_fullTS_PFTs_5 import timedays, outarray, ms

from runs.modelrun_aggTS_RegimeComp import ms1, ms2

from importlib import reload

# TODO:
#  - overlay plots
#  1. show raw forcing as scatter with two separate colors
#  2. show raw forcing with alpha and mean in bold line
#  3. show boxplot next to raw forcing showing model output regimes vs env data
#  4. create a similar plot for model phytoplankton output (bulk)
#  5. create same plot for complex model with PFT differentiation vs HPLC PFT data

params = {'font.size': 7}
matplotlib.rcParams.update(params)

outarray_reg1, timedays_reg1, _ = ms.physics.forcing.NOX.returnModelOut_Regime(outarray,timedays, regime=1)
outarray_reg2, timedays_reg2, _ = ms.physics.forcing.NOX.returnModelOut_Regime(outarray,timedays, regime=2)

timedays_reg1x = np.mod(timedays_reg1, 365)
timedays_reg2x = np.mod(timedays_reg2, 365)

if ms.physics.forcing.NOX.extendstart:
    print(len(outarray), len(timedays))
    outarray = ms.physics.forcing.NOX.returnModelOut_nospinup(outarray)
    timedays = ms.physics.forcing.NOX.returnTimeDays_nospinup(timedays)

nn = ms.nutrients.num
pn = ms.phytoplankton.num
zn = ms.zooplankton.num
dn = ms.detritus.num
nuts = slice(0, nn)
phyto = slice(nn, nn+pn)
zoo = slice(nn + pn, nn + pn + zn)
det = slice(nn + pn + zn, nn + pn + zn + dn)
flux = slice(nn + pn + zn + dn, None)

outindex = nn + pn + zn + dn
print('outindex', outindex)


def plotregimex(ax1, ax1s, ax2, ax3, ax4, ms1, timedays_reg1, outarray_reg1, regime):
    # color vectors
    #colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', 'grey']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [2, 2.5, 4, 5.5]

    ax1.set_title(regime + '- model out')

    dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    dpm_cumsum = np.cumsum(dayspermonth) - np.array(dayspermonth) / 2  # - 15
    # print(timedays_ly)

    # Figure 1
    N = ms1.physics.forcing.verif.fullpadNO2NO3

    Si = ms1.physics.forcing.verif.fullpadSiOH
    # Nmean = ms.physics.forcing.verif.returnMeanVerifPerMonth(N,'NO3_NO2_USF_Box')

    # print(N, np.arange(len(N)))
    # print(N)
    # N
    N_Max = 10  # np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) + np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) * 0.1
    ax1.scatter(N['yday'], N['NO3_NO2_USF_Box'], label='data', alpha=0.5)
    # ax1[zuplot].plot(np.arange(len(N)), N['NO3_NO2_USF_Box'], label='mean data', alpha=1)
    ax1.scatter(timedays_reg1, outarray_reg1[:, nuts][:, 0], c=colors[1], lw=lws[0], alpha=0.6, label='Model', s=0.1)
    ax1.set_ylabel('Nitrate + Nitrite \n' '[µM N]', multialignment='center')
    ax1.set_ylim(0, N_Max)
    ax1.legend(fontsize='x-small')

    Si_Max = 10  # np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) + np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) * 0.1
    ax1s.scatter(N['yday'], Si['SiO4_USF_Box'], c=colors[3], label='data', alpha=0.5)
    # ax1[zuplot].plot(np.arange(len(N)), N['NO3_NO2_USF_Box'], label='mean data', alpha=1)
    ax1s.scatter(timedays_reg1, outarray_reg1[:, nuts][:, 1], c=colors[1], lw=lws[0], alpha=0.6, label='Model', s=0.1)
    ax1s.set_ylabel('Silicate \n' '[µM N]', multialignment='center')
    ax1s.set_ylim(0, Si_Max)
    ax1s.legend(fontsize='x-small')

    ChlConv = True
    # Phyto
    CtoChla = 50  # g/g
    MolarMassC = 12.0107
    CtoNratioPhyto = 6.625
    muMolartoChlaconvfactor = CtoChla / MolarMassC / CtoNratioPhyto  # Chla as mg/m-3 to

    chla = ms1.physics.forcing.verif.fullpadFluorChla

    chla2 = ms1.physics.forcing.verif.fullpadHPLC
    # print(chla)
    ax2.scatter(chla['yday'], chla['Chlorophyll_Box'] * muMolartoChlaconvfactor, label='Fluor Chla', alpha=0.5)
    ax2.scatter(chla2['yday'], chla2['Tchla'] * muMolartoChlaconvfactor / 100, label='HPLC Chla', alpha=0.5)

    chlamean = ms1.physics.forcing.verif.returnMeanVerifPerMonth(chla, 'Chlorophyll_Box')

    chla2mean = ms1.physics.forcing.verif.returnMeanVerifPerMonth(chla2, 'Tchla')

    # ax2[zuplot].plot(dpm_cumsum, np.array(chlamean) * muMolartoChlaconvfactor, label='mean FlChla', alpha=1)
    # ax2[zuplot].plot(dpm_cumsum, np.array(chla2mean) * muMolartoChlaconvfactor / 100, label='mean HPLC', alpha=1)
    Pall = outarray_reg1[:, phyto]
    P_Max = 2  # np.max(Pall) + 0.9 * np.max(Pall)
    # print(Pall)
    for i in range(pn):
        ax2.scatter(timedays_reg1, Pall[:, i], c=colors[i + 3], lw=lws[1], label='PFT ' + str(i + 1), alpha=0.6, s=0.1)

    ax2.scatter(timedays_reg1, np.sum(Pall, axis=1), c=colors[0], lw=1., label='sum(P)', alpha=0.9, s=0.1)
    # ax2.plot(timedays_reg1, np.sum(Pall, axis=1), ':', c=colors[0], lw=1., label='sum(P)', alpha=0.9)

    PFT_all = np.sum(Pall, axis=1)

    PFT_1 = Pall[:, 0]

    PFT_2 = Pall[:, 1]

    '''
    print(np.mean(timedays_ly[100:2400]/365))
    print('PFT1',np.mean(PFT_1[100:2400]), np.var(PFT_1[100:2400]))
    print('PFT2',np.mean(PFT_2[100:2400]), np.var(PFT_2[100:2400]))
    print('PFTsum',np.mean(PFT_all[100:2400]), np.var(PFT_all[100:2400]))

    print(np.mean(timedays_ly[11*365:15*365]/365))
    print('PFT1',np.mean(PFT_1[11*365:15*365]), np.var(PFT_1[11*365:15*365]))
    print('PFT2',np.mean(PFT_2[11*365:15*365]), np.var(PFT_2[11*365:15*365]))
    print('PFTsum',np.mean(PFT_all[11*365:15*365]), np.var(PFT_all[11*365:15*365]))
    '''

    ax2.legend(fontsize='x-small')
    ax2.set_ylabel('Phytoplankton \n' '[µM N]', multialignment='center')
    ax2.set_ylim(0, P_Max)
    # ax2.set_yscale('log')

    # mg dry weight per cubic meter to µM of N
    mggramstograms = 1 / 1000
    Cperdryweight = 0.32
    # Wiebe et al. 1975 : Carbon was 31-33% ofdryweight
    molarmassCarbon = 12.01  # grams per mole
    CtonNratioZoo = 5.625
    mgDWtomuMolarZOO = mggramstograms / Cperdryweight / molarmassCarbon / CtonNratioZoo * 1000  # µM

    ZBM = ms1.physics.forcing.verif.fullpadZoo
    # print(ZBM)
    # ZBMmean = ms.physics.forcing.verif.returnMeanVerifPerMonth(ZBM, 'BIOMASS_200')

    # print(ZBMmean)
    ax3.scatter(ZBM['yday'], ZBM['BIOMASS_200'] * mgDWtomuMolarZOO, label='200 µM', alpha=0.5)
    ax3.scatter(ZBM['yday'], ZBM['BIOMASS_500'] * mgDWtomuMolarZOO,
                label='500 µM', alpha=0.5)

    # ax3[zuplot].plot(dpm_cumsum, np.array(ZBMmean) * muMolartoChlaconvfactor * mgDWtomuMolarZOO, label='mean data',
    #                    alpha=1)

    # Z
    Zall = outarray_reg1[:, zoo]
    Z_Max = 1.5  # np.max(Zall) + 0.1 * np.max(Zall)

    ax3.scatter(timedays_reg1, Zall, c=colors[4], lw=lws[1], label='model', alpha=0.6, s=0.1)
    ax3.set_ylabel('Zooplankton \n' '[µM N]', multialignment='center')
    ax3.tick_params('y')
    ax3.set_ylim(0, Z_Max)
    ax3.legend(fontsize='x-small')
    # ax4[i_plot].set_title('Zooplankton')

    # convert PN in µg/L to µM of Detritus!
    molarmassNitrogen = 14.0067
    mugperlitertomuMolarPN = 1 / molarmassNitrogen  # g/L -> mol/L -> µM

    # TODO: add other verification variables! ensure proper converison..
    PN = ms1.physics.forcing.verif.fullpadPN

    # PNmean = ms.physics.forcing.verif.returnMeanVerifPerMonth(PN, 'PON_ug_kg_Box')

    # print(np.array(PNmean),dpm_cumsum)
    # ax4[zuplot].scatter(PN['yday'], PN['PON_ug_kg_Box'] /14.0067 , label='PN data', alpha=0.5)

    D_Max = 2.5  # np.max(outarray_ly[:, 3]) + 0.2 * np.max(outarray_ly[:, 3])
    # D
    ax4.scatter(timedays_reg1, outarray_reg1[:, det], c=colors[1], lw=lws[0], alpha=0.6, label='model', s=0.1)
    ax4.scatter(PN['yday'], PN['PON_ug_kg_Box'] / 14.0067, label='PON data', alpha=0.5)
    ax4.set_ylabel('Detritus \n' '[µM N]', multialignment='center')
    ax4.set_ylim(0, D_Max)
    ax4.set_xlabel('Day in year')
    ax4.legend(fontsize='x-small')
    # Legend


def plotregimefx(ax1, ax1s, ax2, ax3, ax4, msXX, timedays_regXX, outarray_regXX, regime, timedays_reg1XXxY):

    ax1.set_title(regime + '- model out')

    dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    dpm_cumsum = np.cumsum(dayspermonth) - np.array(dayspermonth) / 2  # - 15
    # print(timedays_ly)


    # color vectors
    #colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', 'grey']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [2, 2.5, 4, 5.5]

    P_Max = 2  # np.max(Pall) + 0.9 * np.max(Pall)
    N_Max = 20  # np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) + np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) * 0.1
    D_Max = 1.5  # np.max(outarray_ly[:, 3]) + 0.2 * np.max(outarray_ly[:, 3])
    Z_Max = 0.5  # np.max(Zall) + 0.1 * np.max(Zall)
    P_Max = 2  # np.max(Pall) + 0.9 * np.max(Pall)

    ax1.set_title(regime + 'forcing')

    ax4.set_xlabel('Day in year')

    NOX = msXX.physics.forcing.NOX.return_interpvalattime(timedays_regXX)
    print(len(NOX))
    #print(outarray_regXX, len(outarray_regXX))
    print(timedays_regXX, len(timedays_regXX))
    #print(timedays_reg1XXxY, len(timedays_reg1XXxY))
    NOXdat = msXX.physics.forcing.NOX.forcingfile

    NN = msXX.physics.forcing.NOX.rawforcing

    # ax1[muplot].scatter(NN['yday'], NN['NO3_NO2_USF_AtDepth'] , label='data', alpha=0.5)

    # print(NOX)
    # print(NOXdat)
    ax1.scatter(timedays_reg1XXxY+1, NOX, c=colors[5], lw=lws[0], alpha=alphas[0], s=0.1)
    # ax1[muplot].scatter(dpm_cumsum, NOXdat[0:12], c=colors[5])
    ax1.set_ylabel('$N_0$ \n' '[µM]', multialignment='center')
    ax1.set_ylim(0., N_Max)
    # ax1[muplot].invert_yaxis()

    SiOX = msXX.physics.forcing.SiOH.return_interpvalattime(timedays_regXX)

    SiOXdat = msXX.physics.forcing.SiOH.forcingfile

    SiX = msXX.physics.forcing.SiOH.rawforcing

    # ax1[muplot].scatter(NN['yday'], NN['NO3_NO2_USF_AtDepth'] , label='data', alpha=0.5)

    # print(NOX)
    # print(NOXdat)
    ax1s.scatter(timedays_reg1XXxY+1, SiOX, c=colors[5], lw=lws[0], alpha=alphas[0], s=0.1)
    # ax1[muplot].scatter(dpm_cumsum, NOXdat[0:12], c=colors[5])
    ax1s.set_ylabel('$Si_0$ \n' '[µM]', multialignment='center')
    ax1s.set_ylim(0., N_Max)
    # ax1[muplot].invert_yaxis()

    MLDRAW = msXX.physics.forcing.X258.rawforcing

    # ax2[muplot].scatter(MLDRAW['yday'], MLDRAW['x258depth'], label='data', alpha=0.5)

    MLD = msXX.physics.forcing.X258.return_interpvalattime(timedays_regXX)
    # MLDderiv = ms.physics.forcing.X258.return_derivattime(timedays_ly)
    MLDdat = msXX.physics.forcing.X258.forcingfile
    MLD_max = 150  # np.max(MLDdat) + 0.1 * np.max(MLDdat)

    # MLD = MLDderiv

    ax2.scatter(timedays_reg1XXxY+1, MLD, c=colors[5], lw=lws[0], alpha=alphas[0], s=0.1)
    # ax2[muplot].scatter(dpm_cumsum, MLDdat[0:12], c=colors[5])
    ax2.set_ylabel('$25.8_{iso}$ \n' '[m]', multialignment='center')
    ax2.set_ylim(0, MLD_max) # 400 for biotrans, 100 for Papa
    ax2.invert_yaxis()

    PARAW = msXX.physics.forcing.PAR.rawforcing
    # print(PARAW)
    # ax3[muplot].scatter(PARAW['yday'], PARAW['PAR'], label='data', alpha=0.5)

    PAR = msXX.physics.forcing.PAR.return_interpvalattime(timedays_regXX)
    PARdat = msXX.physics.forcing.PAR.forcingfile
    PAR_max = 70  # np.max(PAR) + 0.1 * np.max(PAR)
    ax3.scatter(timedays_reg1XXxY+1, PAR, c=colors[5], lw=lws[0], alpha=alphas[0], s=0.1)
    # ax3[muplot].scatter(dpm_cumsum, PARdat[0:12], c=colors[5])
    ax3.set_ylabel('PAR \n' '[E $m^{−2}$ $s^{−1}$]', multialignment='center')
    ax3.set_ylim(0, PAR_max)
    # ax1[muplot].invert_yaxis()

    TmldRAW = msXX.physics.forcing.SST.rawforcing
    # print(TmldRAW)
    # ax4[muplot].scatter(TmldRAW['yday'], TmldRAW['Temperature_Box'], label='data', alpha=0.5)

    Tmld = msXX.physics.forcing.SST.return_interpvalattime(timedays_regXX)
    Tmlddat = msXX.physics.forcing.SST.forcingfile
    Tmld_max = 30  # np.max(Tmld) + 0.1 * np.max(Tmld)
    ax4.scatter(timedays_reg1XXxY+1, Tmld, c=colors[5], lw=lws[0], alpha=alphas[0], s=0.1)
    # ax4[muplot].scatter(dpm_cumsum, Tmlddat[0:12], c=colors[5])
    ax4.set_ylabel('$T_{MLD}$ \n' '[°C]', multialignment='center')
    ax4.set_ylim(18, Tmld_max)
    # ax1[muplot].invert_yaxis()














def plotstuff(ms,outarray, zuplot, muplot, regime, plot):
    numcols = 2
    f1, ([ax1,ex1], [ax1s,ex1s], [ax2,ex2], [ax3,ex3], [ax4,ex4]) = plt.subplots(5, numcols, sharex='col',
                                            gridspec_kw={'height_ratios': [1, 1 , 3, 1, 1]})  # , sharey='row')

    # Defining custom 'xlim' and 'ylim' values.
    xlim = (0, 365)

    plt.setp(([ax1,ex1], [ax1s,ex1s], [ax2,ex2], [ax3,ex3], [ax4,ex4]), xticks=[1, 60, 120, 180, 240, 300, 365], xlim=xlim)

    from matplotlib.ticker import MaxNLocator
    for axe in ([ax1,ex1], [ax1s,ex1s], [ax2,ex2], [ax3,ex3], [ax4,ex4]):
        for i in range(numcols):
            axe[i].get_yaxis().set_major_locator(MaxNLocator(nbins=4))
            axe[i].tick_params(top=True, right=True)

    # PLOTTING
    timedays_ly = timedays
    outarray_ly = outarray

    # color vectors
    #colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', 'grey']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [2, 2.5, 4, 5.5]


    if plot == 'model':
        plotregimex(ax1,ax1s,ax2,ax3,ax4,ms1, timedays_reg1x, outarray_reg1,regime)

        plotregimex(ex1,ex1s,ex2,ex3,ex4,ms2,timedays_reg2x, outarray_reg2,'Regime 2')

    ## PHYSICS ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    elif plot=='physics':
        plotregimefx(ax1,ax1s,ax2,ax3,ax4,ms, timedays_reg1, outarray_reg1,regime, timedays_reg1x)

        plotregimefx(ex1,ex1s,ex2,ex3,ex4,ms,timedays_reg2, outarray_reg2,'Regime 2', timedays_reg2x,)

    plt.show()



plotstuff(ms, outarray, 0, 1, 'Regime 1', plot='model')


plotstuff(ms, outarray, 0, 1, 'Regime 1', plot='physics')

# Defining custom 'xlim' and 'ylim' values.
xlim = (0, 365)

# Setting the values for all axes.

#f1.align_ylabels()


#plt.subplots_adjust(vspace=0, hspace=0)

#plt.tight_layout(rect=[1, 1, 1, 1])
#plt.show()