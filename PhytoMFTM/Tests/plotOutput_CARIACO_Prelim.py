#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time
import pandas

# Fitting
from lmfit import minimize, Parameters, Parameter, report_fit

from Tests.runModel_CARIACO_Prelim import out4P2Z,timedays_model


# make all plots larger and more visible on dark background:
plt.rcParams['figure.figsize'] = [16, 10]
plt.rc_context({'axes.edgecolor':'black', 'xtick.color':'black', 'ytick.color':'black', 'figure.facecolor':'white'})

plt.rcParams['figure.dpi']= 300


# read yearly data (for comparison to model) from Cariaco
NO3NO2 = pandas.read_csv('ValidationData/NO2NO3_above_R1.csv')
SiOH = pandas.read_csv('ValidationData/SiOH_above_R1.csv')
DIATOM = pandas.read_csv('ValidationData/ DIATOMS _above_R1.csv')
COCCO = pandas.read_csv('ValidationData/ COCCOLITHOPHORIDS _above_R1.csv')
DINO = pandas.read_csv('ValidationData/ DINOFLAGELLATES _above_R1.csv')
NANO = pandas.read_csv('ValidationData/ NANOFLAGELLATES _above_R1.csv')
Zoo200BM = pandas.read_csv('ValidationData/200BIOMASS_above_R2.csv')
Zoo500BM = pandas.read_csv('ValidationData/500BIOMASS_above_R2.csv')
# zooplankton biomass is in [mg/m^3 dry weight]
# for now, a 10% N content of the dry weight of all Zooplankton is assumed,
# [Gorsky et al. 1988] "C and N composition of some northwestern Mediterranean zooplankton and micronekton species"
# to convert the values to µM  N  ::  1 mg N/m^3 = 0.071394 μM N

###########--monthly medians---######
NO3NO2 = NO3NO2.assign(month = pandas.to_datetime(NO3NO2['yday'], format='%j').dt.month)
SiOH = SiOH.assign(month = pandas.to_datetime(SiOH['yday'], format='%j').dt.month)
DIATOM = DIATOM.assign(month = pandas.to_datetime(DIATOM['yday'], format='%j').dt.month)
COCCO = COCCO.assign(month = pandas.to_datetime(COCCO['yday'], format='%j').dt.month)
DINO = DINO.assign(month = pandas.to_datetime(DINO['yday'], format='%j').dt.month)
NANO = NANO.assign(month = pandas.to_datetime(NANO['yday'], format='%j').dt.month)
Zoo200BM = Zoo200BM.assign(month = pandas.to_datetime(Zoo200BM['yday'], format='%j').dt.month)
Zoo500BM = Zoo500BM.assign(month = pandas.to_datetime(Zoo500BM['yday'], format='%j').dt.month)

NO3NO2_monthly_median = NO3NO2.groupby('month').median()
SiOH_monthly_median = SiOH.groupby('month').median()
DIATOM_monthly_median = DIATOM.groupby('month').median()
COCCO_monthly_median = COCCO.groupby('month').median()
DINO_monthly_median = DINO.groupby('month').median()
NANO_monthly_median = NANO.groupby('month').median()
Zoo200BM_monthly_median = Zoo200BM.groupby('month').median()
Zoo500BM_monthly_median = Zoo500BM.groupby('month').median()

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

colormap = pandas.DataFrame({"name" : ["MLD", "Si", "N", "Phyto", "MikroZ", "MesoZ", "D", "SST", "PAR"],
                             "color" : ["#1D71B8", "#FFDE00", "#575756", "#009640", "#EA5B0C", "#BE1622", "#B17F4A", "#E94E1B", "#F9B233"]})

colormap["name"] = colormap["name"].apply(lambda x: x.lower())
c = dict(zip(*colormap.values.T))
mcolors.get_named_colors_mapping().update(c)


def plotDATAvsYEARoutput(outarray, pfn, zn, i_plot, title):
    f1, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9) = plt.subplots(9, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    # PLOTTING
    timedays = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    i_plot =0
    # color vectors
    colors = ['#808080' ,'#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    ax1.set_title(title)

    # Figure 1
    # N
    ax1.plot(timedays, outarray_ly[:, 0], c="N", lw=lws[0], alpha=alphas[0], label='Model')
    # N Data
    ax1.scatter(NO3NO2['yday'].values, NO3NO2['NO2NO3'].values, c="N", s=4.3, label='Data')
    if i_plot == 0:
        ax1.set_ylabel('Nitrate \n' '[µM]', multialignment='center', fontsize=10)
    # ax1[i_plot].set_ylim(-0.1, 5)

    # Si
    ax2.plot(timedays, outarray_ly[:, 1], c="Si", lw=lws[0], alpha=alphas[0])
    # Si Data
    ax2.scatter(SiOH['yday'].values, SiOH['SiOH'].values, c="Si", s=4.3)
    if i_plot == 0:
        ax2.set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    # ax2[i_plot].set_ylim(-0.1, 12)

    # Phyto i=2
    # ax3.plot(timedays, sum([outarray_ly[:, 3 + zn + i] for i in range(pfn)]), c=colors[4], lw=lws[1])
    ax3.plot(timedays, outarray_ly[:, 3 + zn + 0], c="Phyto")
    #ax3_tx = ax3.twinx()
    ax3.scatter(DIATOM['yday'].values, DIATOM['abundance'].values * 2 / 1000, c="Phyto", s=4.3)
    # 2. pmol per cell [Marchetti and Harrison 2007]
    # values of abundance are organisms per ml
    # to convert pmol per cell per ml to µM = * abundance / 1000

    ax4.plot(timedays, outarray_ly[:, 3 + zn + 1], c="Phyto")
    #ax4_tx = ax4.twinx()
    ax4.scatter(COCCO['yday'].values, COCCO['abundance'].values * 1 / 1000, c="Phyto", s=4.3)
    # 1 pmol per cell [Aksnes et al. 1994]

    ax5.plot(timedays, outarray_ly[:, 3 + zn + 2], c="Phyto")
    #ax5_tx = ax5.twinx()
    ax5.scatter(DINO['yday'].values, DINO['abundance'].values * 1.5 / 1000, c="Phyto", s=4.3)
    # 1.5 pmol per cell [Li et al. 2016]

    ax6.plot(timedays, outarray_ly[:, 3 + zn + 3], c="Phyto")
    #ax6_tx = ax6.twinx()
    ax6.scatter(NANO['yday'].values, NANO['abundance'].values * 0.1 / 1000, c="Phyto", s=4.3)
    # 30 fmol per cell > 0.03 pmol per cell [Maat et al. 2014]

    #ax3_tx.set_ylim(-0.1, 0.8)
    #ax4_tx.set_ylim(-0.1, 0.8)
    #ax5_tx.set_ylim(-0.1, 0.8)
    #ax6_tx.set_ylim(-0.1, 0.8)

    if i_plot == 0:
        ax3.set_ylabel('Diatom \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax4.set_ylabel('Coccs \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax5.set_ylabel('Dinos \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax6.set_ylabel('Nano \n' '[µM N]', multialignment='center', fontsize=10)
    #ax3.set_ylim(-0.1, 0.8)
    #ax4.set_ylim(-0.1, 0.8)
    #ax5.set_ylim(-0.1, 0.8)
    #ax6.set_ylim(-0.1, 0.8)

    # ax3[i_plot].set_title('Phy Biomass & ChlA Data')

    # Z
    i=3
    # ax4.plot(timedays, sum([outarray_ly[:, 3 + i] for i in range(zn)]), c=colors[4], lw=lws[1])
    # [ax4.plot(timedays, outarray_ly[:, 3 + i], c=colors[i + 1], lw=lws[0], alpha=alphas[0]) for i in range(zn)]
    ax7.plot(timedays, outarray_ly[:, 3 + 0], c="MikroZ", lw=lws[0], alpha=alphas[0])
    ax7.scatter(Zoo200BM['yday'].values, Zoo200BM['abundance'].values *0.1*0.071394, c="MikroZ", s=4.3)

    ax8.plot(timedays, outarray_ly[:, 3 + 1], c="MesoZ", lw=lws[0], alpha=alphas[0])
    ax8.scatter(Zoo500BM['yday'].values, Zoo500BM['abundance'].values *0.1*0.071394, c="MesoZ", s=4.3)
    if i_plot == 0:
        ax7.set_ylabel('Mikro Z \n' '[µM N]', multialignment='center', fontsize=9)
    if i_plot == 0:
        ax8.set_ylabel('Meso Z \n' '[µM N]', multialignment='center', fontsize=9)
    ax7.tick_params('y', labelsize=10)

    # ax4[i_plot].set_title('Zooplankton')
    # ax4[i_plot].set_ylim(0, 0.62)

    # D
    ax9.plot(timedays, outarray_ly[:, 2], c="D", lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax9.set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    # ax5[i_plot].set_title('Detritus')
    # ax5[i_plot].set_ylim(0,0.15)

    ax9.set_xlabel('Day in year')
    # Legend
    f1.align_ylabels()
    plt.margins(x=0)
    #adjustFigAspect(fig, aspect=.5)
    plt.tight_layout()
    plt.savefig('FirstNaiveOutputCARIACO.png')





plotDATAvsYEARoutput(out4P2Z, 4, 2, 1, '4P2Z')