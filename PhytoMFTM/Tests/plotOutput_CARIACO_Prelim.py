#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time
import pandas

# Fitting
from lmfit import minimize, Parameters, Parameter, report_fit

from Tests.runModel_CARIACO_Prelim import timedays_model, out5P2Zconstant # , out5P2Z,


# make all plots larger and more visible on dark background:
#plt.rcParams['figure.figsize'] = [16, 10]
#plt.rc_context({'axes.edgecolor':'black', 'xtick.color':'black', 'ytick.color':'black', 'figure.facecolor':'white'})

#plt.rcParams['figure.dpi']= 300


# read yearly data (for comparison to model) from Cariaco
NO3NO2 = pandas.read_csv('ValidationData/Old/NewVerifDATA/NO2NO3_r1.csv')
SiOH = pandas.read_csv('ValidationData/Old/NewVerifDATA/SiOH_r1.csv')

DIATOM = pandas.read_csv('ValidationData/Old/NewVerifDATA/Diatoms_r1.csv')
HAPTO = pandas.read_csv('ValidationData/Old/NewVerifDATA/Hapto_r1.csv')
DINO = pandas.read_csv('ValidationData/Old/NewVerifDATA/Dino_r1.csv')
CYANO = pandas.read_csv('ValidationData/Old/NewVerifDATA/Cyano_r1.csv')
OTHERS = pandas.read_csv('ValidationData/Old/NewVerifDATA/Others_r1.csv')

Zoo200BM = pandas.read_csv('ValidationData/Old/200BIOMASS_above_R2.csv')
Zoo500BM = pandas.read_csv('ValidationData/Old/500BIOMASS_above_R2.csv')

PN = pandas.read_csv('ValidationData/Old/NewVerifDATA/PN_r1.csv')
# zooplankton biomass is in [mg/m^3 dry weight]
# for now, a 10% N content of the dry weight of all Zooplankton is assumed,
# [Gorsky et al. 1988] "C and N composition of some northwestern Mediterranean zooplankton and micronekton species"
# to convert the values to µM  N  ::  1 mg N/m^3 = 0.071394 μM N

###########--monthly medians---######
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
#DIATOM_monthly_median = DIATOM.groupby('month').median()
#COCCO_monthly_median = COCCO.groupby('month').median()
#DINO_monthly_median = DINO.groupby('month').median()
#NANO_monthly_median = NANO.groupby('month').median()
Zoo200BM_monthly_median = Zoo200BM.groupby('month').median()
Zoo500BM_monthly_median = Zoo500BM.groupby('month').median()

# COLOR SCHEMwE

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
    i_plot =0
    # color vectors
    colors = ['#808080', '#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

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

    ax1[0].plot(timedays, outarray_ly[:, 3 + zn + 0], c="Phyto")
    ax1[0].set_ylabel('Diatoms \n' '[µM N]', multialignment='center', fontsize=9)
    ax1_tx = ax1[0].twinx()
    ax1_tx.scatter(DIATOM['yday'].values, DIATOM['Value'].values / 100, c="Phyto", s=4.3)
    # 2. pmol per cell [Marchetti and Harrison 2007]
    # values of abundance are organisms per ml
    # to convert pmol per cell per ml to µM = * abundance / 1000

    ax2[0].plot(timedays, outarray_ly[:, 3 + zn + 1], c="Phyto")
    ax2[0].set_ylabel('Hapto \n' '[µM N]', multialignment='center', fontsize=9)
    ax2_tx = ax2[0].twinx()
    ax2_tx.scatter(HAPTO['yday'].values, HAPTO['Value'].values / 100, c="Phyto", s=4.3)

    ax3[0].plot(timedays, outarray_ly[:, 3 + zn + 2], c="Phyto")
    ax3[0].set_ylabel('Cyano \n' '[µM N]', multialignment='center', fontsize=9)
    ax3_tx = ax3[0].twinx()
    ax3_tx.scatter(CYANO['yday'].values, CYANO['Value'].values / 100, c="Phyto", s=4.3)

    ax4[0].plot(timedays, outarray_ly[:, 3 + zn + 3], c="Phyto")
    ax4[0].set_ylabel('Dino \n' '[µM N]', multialignment='center', fontsize=9)
    ax4_tx = ax4[0].twinx()
    ax4_tx.scatter(DINO['yday'].values, DINO['Value'].values / 100, c="Phyto", s=4.3)

    ax5[0].plot(timedays, outarray_ly[:, 3 + zn + 4], c="Phyto")
    ax5[0].set_ylabel('Others \n' '[µM N]', multialignment='center', fontsize=9)
    ax5_tx = ax5[0].twinx()
    ax5_tx.scatter(OTHERS['yday'].values, OTHERS['Value'].values / 100, c="Phyto", s=4.3)

    i=3
    ax3[1].plot(timedays, outarray_ly[:, 3 + 0], c="MikroZ", lw=lws[0], alpha=alphas[0])
    ax3[1].scatter(Zoo200BM['yday'].values, Zoo200BM['Value'].values *0.1*0.071394, c="MikroZ", s=4.3)

    ax4[1].plot(timedays, outarray_ly[:, 3 + 1], c="MesoZ", lw=lws[0], alpha=alphas[0])
    ax4[1].scatter(Zoo500BM['yday'].values, Zoo500BM['Value'].values *0.1*0.071394, c="MesoZ", s=4.3)

    ax3[1].set_ylabel('Mikro Z \n' '[µM N]', multialignment='center', fontsize=9)

    ax4[1].set_ylabel('Meso Z \n' '[µM N]', multialignment='center', fontsize=9)

    ax2[1].tick_params('y', labelsize=10)

    # D
    ax5[1].plot(timedays, outarray_ly[:, 2], c="De", lw=lws[0], alpha=alphas[0])
    ax5[1].set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)
    ax51_tx = ax5[1].twinx()
    ax51_tx.scatter(PN['yday'].values, PN['Value'].values, c="De", s=4.3)

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
    plt.margins(x=0)
    #adjustFigAspect(fig, aspect=.5)
    #plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    #plt.savefig('FirstNaiveOutputCARIACO.png')


def plotMODELFORCINGoutput(outarray, pfn, zn):
    f2, (ax1, ax2, ax22, ax3, ax4) = plt.subplots(5, 1, sharex='col')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    # PLOTTING
    timedays = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    i_plot =0
    # color vectors
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # Figure 1
    # N
    ax1.plot(timedays, outarray_ly[:, 3+zn+pfn], c="Ni", lw=lws[0], alpha=alphas[0], label='Model')
    # N Data
    ax1.set_ylim(0, 101)
    ax1.invert_yaxis()
    ax1.set_ylabel('MLD \n' '[m]', multialignment='center', fontsize=10)

    # Si
    ax2.plot(timedays, outarray_ly[:, 3+zn+pfn+2], c="Ni", lw=lws[0], alpha=alphas[0], label='K')
    ax22.plot(timedays, outarray_ly[:, 3+zn+pfn+4], c="MLD", lw=lws[0], alpha=alphas[0], label='U')
    # Si Data
    ax2.set_ylabel('K \n' '[$d^{-1}$]', multialignment='center', fontsize=10)
    ax22.set_ylabel('U \n' '[$d^{-1}$]', multialignment='center', fontsize=10)
    ax2.legend()

    ax3.plot(timedays, outarray_ly[:, 3+zn+pfn+1], c="Ni", lw=lws[0], alpha=alphas[0])
    ax3.set_ylabel('NMixing \n' '[µM $d^{-1}$]', multialignment='center', fontsize=10)

    ax4.plot(timedays, outarray_ly[:, 3+zn+pfn+3], c="Phyto", lw=lws[0], alpha=alphas[0])
    ax4.set_ylabel('Gains \n' '[µM $d^{-1}$]', multialignment='center', fontsize=10)
    ax4.set_xlabel('Day in year')
    # Legend
    f2.align_ylabels()
    #f1.delaxes(ax = ax1[0])
    #plt.margins(x=0)
    #adjustFigAspect(fig, aspect=.5)
    #plt.tight_layout()

    plt.subplots_adjust(hspace=0,wspace=0)
    #plt.savefig('FirstNaiveOutputCARIACO.png')


#plotDATAvsYEARoutput(out5P2Z, 5, 2, 1, '5P2Z')

#plotMODELFORCINGoutput(out5P2Z, 5, 2)

plotDATAvsYEARoutput(out5P2Zconstant, 5, 2, 1, '5P2Z')

plotMODELFORCINGoutput(out5P2Zconstant, 5, 2)
