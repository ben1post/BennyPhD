#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.integrate import odeint
import time
import pandas

# Fitting
from lmfit import minimize, Parameters, Parameter, report_fit

from PhDPropPlots.runModel_BioticEnvMixing import outarray, timedays_model #, out5P2Zconstant, out5P2Z, out5P2Z_2


out5P2Z_3 = outarray
# make all plots larger and more visible on dark background:
#plt.rcParams['figure.figsize'] = [16, 10]
#plt.rc_context({'axes.edgecolor':'black', 'xtick.color':'black', 'ytick.color':'black', 'figure.facecolor':'white'})

#plt.rcParams['figure.dpi']= 300


# read yearly data (for comparison to model) from Cariaco
NO3NO2 = pandas.read_csv('ValidationData/NewVerifDATA/NO2NO3_r1.csv')
SiOH = pandas.read_csv('ValidationData/NewVerifDATA/SiOH_r1.csv')

DIATOM = pandas.read_csv('ValidationData/NewVerifDATA/Diatoms_r1.csv')
HAPTO = pandas.read_csv('ValidationData/NewVerifDATA/Hapto_r1.csv')
DINO = pandas.read_csv('ValidationData/NewVerifDATA/Dino_r1.csv')
CYANO = pandas.read_csv('ValidationData/NewVerifDATA/Cyano_r1.csv')
OTHERS = pandas.read_csv('ValidationData/NewVerifDATA/Others_r1.csv')

Zoo200BM = pandas.read_csv('ValidationData/200BIOMASS_above_R2.csv')
Zoo500BM = pandas.read_csv('ValidationData/500BIOMASS_above_R2.csv')

PN = pandas.read_csv('ValidationData/NewVerifDATA/PN_r1.csv')
# zooplankton biomass is in [mg/m^3 dry weight]
# for now, a 10% N content of the dry weight of all Zooplankton is assumed,
# [Gorsky et al. 1988] "C and N composition of some northwestern Mediterranean zooplankton and micronekton species"
# to convert the values to µM  N  ::  1 mg N/m^3 = 0.071394 μM N

NO3NO2 = NO3NO2.assign(month = pandas.to_datetime(NO3NO2['yday'], format='%j').dt.month)
SiOH = SiOH.assign(month = pandas.to_datetime(SiOH['yday'], format='%j').dt.month)

DIATOM = DIATOM.assign(month = pandas.to_datetime(DIATOM['yday'], format='%j').dt.month)
HAPTO = HAPTO.assign(month = pandas.to_datetime(HAPTO['yday'], format='%j').dt.month)
DINO = DINO.assign(month = pandas.to_datetime(DINO['yday'], format='%j').dt.month)
CYANO = CYANO.assign(month = pandas.to_datetime(CYANO['yday'], format='%j').dt.month)
OTHERS = OTHERS.assign(month = pandas.to_datetime(OTHERS['yday'], format='%j').dt.month)

Zoo200BM = Zoo200BM.assign(month = pandas.to_datetime(Zoo200BM['yday'], format='%j').dt.month)
Zoo500BM = Zoo500BM.assign(month = pandas.to_datetime(Zoo500BM['yday'], format='%j').dt.month)

PN = PN.assign(month = pandas.to_datetime(PN['yday'], format='%j').dt.month)


NO3NO2_monthly_median = NO3NO2.groupby('month').median()
SiOH_monthly_median = SiOH.groupby('month').median()

DIATOM_monthly_median = DIATOM.groupby('month').median()
HAPTO_monthly_median = HAPTO.groupby('month').median()
DINO_monthly_median = DINO.groupby('month').median()
CYANO_monthly_median = CYANO.groupby('month').median()
OTHERS_monthly_median = OTHERS.groupby('month').median()

Zoo200BM_monthly_median = Zoo200BM.groupby('month').median()
Zoo500BM_monthly_median = Zoo500BM.groupby('month').median()

PN_monthly_median = PN.groupby('month').median()

# COLOR SCHEME

# color scheme

# MLD
# FFDE00 Silicate
# 575756 Nitrate
# 009640 Phytoplankton
# EA5B0C Mikrozooplankton
# BE1622 Mesozooplankton
# B17F4A Detritus
# DEDC00 SST
# F9B233 PAR

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

colormap = pandas.DataFrame({"name" : ["MLD", "Si", "Ni", "Phyto", "MikroZ", "MesoZ", "De", "SST", "PAR"],
                             "color" : ["#1D71B8", "#FFDE00", "#575756", "#009640", "#EA5B0C", "#BE1622", "#B17F4A", "#E94E1B", "#F9B233"]})

colormap["name"] = colormap["name"].apply(lambda x: x.lower())
c = dict(zip(*colormap.values.T))
mcolors.get_named_colors_mapping().update(c)


def plotDATAvsYEARoutput(outarray, pfn, zn, i_plot, title):
    f1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 2, sharex='col')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    # PLOTTING
    timedays = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    i_plot = 0
    # color vectors
    colors = ['#808080', '#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    modeldepth = 100
    CtoNratioPhyto = 6.625
    CtoChla = 50
    muMolartoChlaconvfactor = modeldepth * CtoNratioPhyto / CtoChla # µM N m^-2 * C N^-1 / C chla^-1 = µM chla m^-2

    # mg dry weight per cubic meter to µM of N
    mggramstograms = 1/1000
    Cperdryweight = 0.32
    # Wiebe et al. 1975 : Carbon was 31-33% ofdryweight
    molarmassCarbon = 12.01 #grams per mole
    CtonNratioZoo = 5.625
    mgDWtomuMolarZOO = mggramstograms / Cperdryweight / molarmassCarbon / CtonNratioZoo * 1000 # µM

    # convert PN in µg/L to µM of Detritus!
    molarmassNitrogen = 14.0067
    mugperlitertomuMolarPN = 1 / molarmassNitrogen # g/L -> mol/L -> µM

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    # Figure 1
    # N
    ax1[1].plot(timedays, outarray_ly[:, 0], c="Ni", lw=lws[0], alpha=alphas[0], label='Model')
    # N Data
    ax1[1].scatter(NO3NO2['yday'].values, NO3NO2['Value'].values, c="Ni", s=4.3, label='Data')
    ax1[1].set_ylabel('Nitrate \n' '[µM]', multialignment='center', fontsize=10)

    # Si
    ax2[1].plot(timedays, outarray_ly[:, 1], c="Si", lw=lws[0], alpha=alphas[0])
    # Si Data
    ax2[1].scatter(SiOH['yday'].values, SiOH['Value'].values, c="Si", s=4.3)

    ax2[1].set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)

    ax1[0].plot(timedays, outarray_ly[:, 3 + zn + 0] * muMolartoChlaconvfactor, c="Phyto")
    ax1[0].set_ylabel('Diatoms \n' '[mg Chl m$^{-2}$]', multialignment='center', fontsize=9)
    #ax1_tx = ax1[0].twinx()
    ax1[0].scatter(DIATOM['yday'].values, DIATOM['Value'].values, c="Phyto", s=4.3)
    # 2. pmol per cell [Marchetti and Harrison 2007]
    # values of abundance are organisms per ml
    # to convert pmol per cell per ml to µM = * abundance / 1000

    ax2[0].plot(timedays, outarray_ly[:, 3 + zn + 1] * muMolartoChlaconvfactor, c="Phyto")
    ax2[0].set_ylabel('Hapto \n' '[mg Chl m$^{-2}$]', multialignment='center', fontsize=9)
    #ax2_tx = ax2[0].twinx()
    ax2[0].scatter(HAPTO['yday'].values, HAPTO['Value'].values, c="Phyto", s=4.3)

    ax3[0].plot(timedays, outarray_ly[:, 3 + zn + 2] * muMolartoChlaconvfactor, c="Phyto")
    ax3[0].set_ylabel('Cyano \n' '[mg Chl m$^{-2}$]', multialignment='center', fontsize=9)
    #ax3_tx = ax3[0].twinx()
    ax3[0].scatter(CYANO['yday'].values, CYANO['Value'].values, c="Phyto", s=4.3)

    ax4[0].plot(timedays, outarray_ly[:, 3 + zn + 3] * muMolartoChlaconvfactor, c="Phyto")
    ax4[0].set_ylabel('Dino \n' '[mg Chl m$^{-2}$]', multialignment='center', fontsize=9)
    #ax4_tx = ax4[0].twinx()
    ax4[0].scatter(DINO['yday'].values, DINO['Value'].values, c="Phyto", s=4.3)

    ax5[0].plot(timedays, outarray_ly[:, 3 + zn + 4] * muMolartoChlaconvfactor, c="Phyto")
    ax5[0].set_ylabel('Others \n' '[mg Chl m$^{-2}$]', multialignment='center', fontsize=9)
    #ax5_tx = ax5[0].twinx()
    ax5[0].scatter(OTHERS['yday'].values, OTHERS['Value'].values, c="Phyto", s=4.3)


    i=3
    ax3[1].plot(timedays, outarray_ly[:, 3 + 0], c="MikroZ", lw=lws[0], alpha=alphas[0])
    ax3[1].scatter(Zoo200BM['yday'].values, Zoo200BM['Value'].values * mgDWtomuMolarZOO, c="MikroZ", s=4.3)

    ax4[1].plot(timedays, outarray_ly[:, 3 + 1], c="MesoZ", lw=lws[0], alpha=alphas[0])
    ax4[1].scatter(Zoo500BM['yday'].values, Zoo500BM['Value'].values * mgDWtomuMolarZOO, c="MesoZ", s=4.3)

    ax3[1].set_ylabel('Mikro Z \n' '[µM N]', multialignment='center', fontsize=9)

    ax4[1].set_ylabel('Meso Z \n' '[µM N]', multialignment='center', fontsize=9)

    ax2[1].tick_params('y', labelsize=10)

    # D
    ax5[1].plot(timedays, outarray_ly[:, 2], c="De", lw=lws[0], alpha=alphas[0])
    ax5[1].set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)
    #ax51_tx = ax5[1].twinx()
    ax5[1].scatter(PN['yday'].values, PN['Value'].values * mugperlitertomuMolarPN, c="De", s=4.3)

    # ax5[i_plot].set_title('Detritus')
    # ax5[i_plot].set_ylim(0,0.15)
    ax1[1].yaxis.set_label_position('right')
    ax1[1].yaxis.tick_right()
    ax2[1].yaxis.set_label_position('right')
    ax2[1].yaxis.tick_right()
    ax3[1].yaxis.set_label_position('right')
    ax3[1].yaxis.tick_right()
    ax4[1].yaxis.set_label_position('right')
    ax4[1].yaxis.tick_right()
    ax5[1].yaxis.set_label_position('right')
    ax5[1].yaxis.tick_right()

    #plt.yticks(np.arange(0, 3, .1))

    ax5[0].set_xlabel('Day in year')
    ax5[1].set_xlabel('Day in year')
    # Legend
    f1.align_ylabels()
    #f1.delaxes(ax = ax1[0])
    #plt.margins(x=0)
    #adjustFigAspect(fig, aspect=.5)
    plt.tight_layout()
    #plt.savefig('FirstNaiveOutputCARIACO.png')


def plotMODELFORCINGoutput(outarray, outarray2, outarray3, pfn, zn):
    f2, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 3, gridspec_kw = {'height_ratios':[1, 3, 1, 3, 1]}, sharex='col')#, sharey='row')


    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    # PLOTTING
    timedays = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    outarray_ly_2 = outarray2[1460:1825]
    outarray_ly_3 = outarray3[1460:1825]

    i_plot = 0
    # color vectors
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    ##############################################################################################################

    # Figure 3
    P_3 = outarray_ly_3[:, 3 + zn + pfn + 7]
    PGains_3 = outarray_ly_3[:, 3 + zn + pfn + 3]
    PGrazed_3 = outarray_ly_3[:, 3 + zn + pfn + 8]
    PMix_3 = outarray_ly_3[:, 3 + zn + pfn + 9]
    PSink_3 = outarray_ly_3[:, 3 + zn + pfn + 10]
    PMort_3 = outarray_ly_3[:, 3 + zn + pfn + 17]

    ax1[2].set_title('Constant $N_0$ + bbox = 100m', fontsize=10)

    ax1[2].plot(timedays, P_3, label='P', color='Phyto')
    ax1[2].set_ylim(0,4)

    ax2[2].stackplot(timedays, PGains_3, labels=['P Gains'], baseline='zero')
    ax2[2].stackplot(timedays, -PGrazed_3, -PMix_3, -PSink_3, -PMort_3, labels=['P Grazed', 'P Mix', 'P Sink', 'P Mortality'], baseline='zero')
    ax2[2].plot(timedays, PGains_3 - PGrazed_3 - PMix_3 - PSink_3 - PMort_3, label='Total Flux', color='black')
    ax2[2].set_ylim(-1., 1.)

    ax2[2].legend(loc='lower right', fontsize=6)

    Z_3 = outarray_ly_3[:, 3 + zn + pfn + 11]
    ZAssimGraz_3 = outarray_ly_3[:, 3 + zn + pfn + 12]
    ZInterGraz_3 = outarray_ly_3[:, 3 + zn + pfn + 13]
    ZMix_3 = outarray_ly_3[:, 3 + zn + pfn + 14]
    ZMort_3 = outarray_ly_3[:, 3 + zn + pfn + 15]
    ZHighOrdPred_3 = outarray_ly_3[:, 3 + zn + pfn + 16]

    ax3[2].plot(timedays, Z_3, label='Z', color='MikroZ')
    ax3[2].set_ylim(0,4)

    ax4[2].stackplot(timedays, ZAssimGraz_3, -ZMix_3, labels=['Z Assim Graz', 'Z Mix'], baseline='zero')
    ax4[2].stackplot(timedays, -ZInterGraz_3, -ZMort_3, -ZHighOrdPred_3, labels=['Z Inter Graz', 'Z Mortality', 'Z HighOrdPred'], baseline='zero')
    ax4[2].plot(timedays, ZAssimGraz_3 - ZInterGraz_3 - ZMix_3 - ZMort_3 - ZHighOrdPred_3, label='Total Flux', color='black')
    ax4[2].set_ylim(-0.5, 0.5)

    ax4[2].legend(loc='lower right', fontsize=6)

    D_3 = outarray_ly_3[:, 2]

    ax5[2].plot(timedays, D_3, label='D', color='De')
    ax5[2].set_ylim(0,20)

    ax5[2].set_xlabel('Time [days]')
    ax5[2].set_xlim(1,365)
    ax5[2].set_xticks([1, 90, 180, 270, 365])
    plt.subplots_adjust(hspace=0)
    # Legend
    f2.align_ylabels()
    ##############################################################################################################

    # Figure 2
    P_2 = outarray_ly_2[:, 3 + zn + pfn + 7]
    PGains_2 = outarray_ly_2[:, 3 + zn + pfn + 3]
    PGrazed_2 = outarray_ly_2[:, 3 + zn + pfn + 8]
    PMix_2 = outarray_ly_2[:, 3 + zn + pfn + 9]
    PSink_2 = outarray_ly_2[:, 3 + zn + pfn + 10]
    PMort_2 = outarray_ly_2[:, 3 + zn + pfn + 17]

    ax1[1].set_title('Constant $N_0$ + bbox = MLD', fontsize=10)

    ax1[1].plot(timedays, P_2, label='P', color='Phyto')
    ax1[1].set_ylim(0, 4)

    ax2[1].stackplot(timedays, PGains_2, labels=['P Gains'], baseline='zero')
    ax2[1].stackplot(timedays, -PGrazed_2, -PMix_2, -PSink_2, -PMort_2, labels=['P Grazed', 'P Mix', 'P Sink', 'P Mortality'], baseline='zero')
    ax2[1].plot(timedays, PGains_2 - PGrazed_2 - PMix_2 - PSink_2 -PMort_2, label='Total Flux', color='black')
    ax2[1].set_ylim(-1., 1.)

    ax2[1].legend(loc='lower right', fontsize=6)

    Z_2 = outarray_ly_2[:, 3 + zn + pfn + 11]
    ZAssimGraz_2 = outarray_ly_2[:, 3 + zn + pfn + 12]
    ZInterGraz_2 = outarray_ly_2[:, 3 + zn + pfn + 13]
    ZMix_2 = outarray_ly_2[:, 3 + zn + pfn + 14]
    ZMort_2 = outarray_ly_2[:, 3 + zn + pfn + 15]
    ZHighOrdPred_2 = outarray_ly_2[:, 3 + zn + pfn + 16]

    ax3[1].plot(timedays, Z_2, label='Z', color='MikroZ')
    ax3[1].set_ylim(0,4)

    ax4[1].stackplot(timedays, ZAssimGraz_2, -ZMix_2, labels=['Z Assim Graz', 'Z Mix'], baseline='zero')
    ax4[1].stackplot(timedays, -ZInterGraz_2, -ZMort_2, -ZHighOrdPred_2, labels=['Z Inter Graz', 'Z Mortality', 'Z HighOrdPred'], baseline='zero')
    ax4[1].plot(timedays, ZAssimGraz_2 - ZInterGraz_2 - ZMix_2 - ZMort_2 - ZHighOrdPred_2, label='Total Flux', color='black')
    ax4[1].set_ylim(-0.5, 0.5)

    ax4[1].legend(loc='lower right', fontsize=6)

    D_2 = outarray_ly_2[:, 2]

    ax5[1].plot(timedays, D_2, label='D', color='De')
    ax5[1].set_ylim(0,20)

    ax5[1].set_xlabel('Time [days]')
    ax5[1].set_xlim(1,365)
    ax5[1].set_xticks([1, 90, 180, 270, 365])
    plt.subplots_adjust(hspace=0)
    # Legend
    f2.align_ylabels()
    ##############################################################################################################

    # Figure 1
    P = outarray_ly[:, 3 + zn + pfn + 7]
    PGains = outarray_ly[:, 3 + zn + pfn + 3]
    PGrazed = outarray_ly[:, 3 + zn + pfn + 8]
    PMix = outarray_ly[:, 3 + zn + pfn + 9]
    PSink = outarray_ly[:, 3 + zn + pfn + 10]
    PMort = outarray_ly[:, 3 + zn + pfn + 17]

    ax1[0].set_title('Variable $N_0$ + bbox = MLD', fontsize=10)

    ax1[0].plot(timedays, P, label='P', color='Phyto')
    ax1[0].set_ylim(bottom=0)

    ax2[0].stackplot(timedays, PGains, labels=['P Gains'], baseline='zero')
    ax2[0].stackplot(timedays, -PGrazed, -PMix, -PSink, -PMort, labels=['P Grazed', 'P Mix', 'P Sink', 'P Mortality'], baseline='zero')
    ax2[0].plot(timedays, PGains - PGrazed - PMix - PSink -PMort, label='Total Flux', color='black')
    ax2[0].set_ylim(-0.1, 0.1)

    ax2[0].legend(loc='lower right', fontsize=6)

    Z = outarray_ly[:, 3 + zn + pfn + 11]
    ZAssimGraz = outarray_ly[:, 3 + zn + pfn + 12]
    ZInterGraz = outarray_ly[:, 3 + zn + pfn + 13]
    ZMix = outarray_ly[:, 3 + zn + pfn + 14]
    ZMort = outarray_ly[:, 3 + zn + pfn + 15]
    ZHighOrdPred = outarray_ly[:, 3 + zn + pfn + 16]

    ax3[0].plot(timedays, Z, label='Z', color='MikroZ')
    ax3[0].set_ylim(0, 1)

    ax4[0].stackplot(timedays, ZAssimGraz, -ZMix, labels=['Z Assim Graz', 'Z Mix'], baseline='zero')
    ax4[0].stackplot(timedays, -ZInterGraz, -ZMort, -ZHighOrdPred, labels=['Z Inter Graz', 'Z Mortality', 'Z HighOrdPred'], baseline='zero')
    ax4[0].plot(timedays, ZAssimGraz - ZInterGraz - ZMix - ZMort - ZHighOrdPred, label='Total Flux', color='black')
    ax4[0].set_ylim(-0.05, 0.05)

    ax4[0].legend(loc='lower right', fontsize=6)

    D = outarray_ly[:, 2]

    ax5[0].plot(timedays, D, label='D', color='De')
    ax5[0].set_ylim(bottom=0)

    ax5[0].set_xlabel('Time [days]')
    ax5[0].set_xlim(1,365)
    ax5[0].set_xticks([1, 90, 180, 270, 365])
    plt.subplots_adjust(hspace=0)
    # Legend
    f2.align_ylabels()


    ax1[0].set_ylabel('Phyto (Sum) \n' '[µM]', multialignment='center', fontsize=10)

    ax2[0].set_ylabel('Phyto Fluxes \n' '[µM $d^{-1}$]', multialignment='center', fontsize=10)

    ax3[0].set_ylabel('Zoo (Sum) \n' '[µM]', multialignment='center', fontsize=10)

    ax4[0].set_ylabel('Zoo Fluxes \n' '[µM $d^{-1}$]', multialignment='center', fontsize=10)

    ax5[0].set_ylabel('Detritus \n' '[µM]', multialignment='center', fontsize=10)

    plt.margins(x=0)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.2)
    #plt.savefig('PhysicalEnv_CARIACO_firstTESTS.pdf', dpi=300)


plotDATAvsYEARoutput(out5P2Z_3, 5, 2, 1, 'hey')

#plotMODELFORCINGoutput(out5P2Z, out5P2Z_2, out5P2Z_3, 5, 2)