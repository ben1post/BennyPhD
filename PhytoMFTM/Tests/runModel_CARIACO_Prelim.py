#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import time
import pandas

# Fitting
from lmfit import minimize, Parameters, Parameter, report_fit

# loading Modules
from PhytoMFTM.ModelClasses import Plankton
import PhytoMFTM.ModelCode as mc


#params.add('xvar', value=0.50, min=0, max=1)
#params.add('yvar', expr='1.0 - xvar')

# set up basic parameters of model:
standardparams = Parameters()

# number of phytoplankton func types
# standardparams.add('pfun_num', value=4, vary=False)
# number of zooplankton groups
# standardparams.add('zoo_num', value=2, vary=False)

# mld - related
standardparams.add('kappa', value=0.1, vary=False) # min=0.09, max=0.11) # Diffusive mixing across thermocline (m*d^-1)
standardparams.add('deltaD_N', value=0.05, vary=False)   # Nitrate Mineralization rate (d^-1)

standardparams.add('kw', value=0.1, vary=False)     # Light attenuation constant (m^-1)

standardparams.add('v', value=0.11, vary=False)      # Sinking of Phytoplankton from Mixed Layer
standardparams.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)

# p - related
standardparams.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)

standardparams.add('ratioSi', value=0, vary=False)  # Silicate ratio ## FALLBACK PARAM FOR OTHER PFTs
standardparams.add('U_N', value=0, vary=False)    # Nitrate Half Saturation Constant
standardparams.add('U_Si', value=0, vary=False)   # Silicate Half Saturation Constant
standardparams.add('muP', value=0, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 1 (e.g. DIATOMS)
ptype1 = Parameters()
ptype1.add('pt1_ratioSi', value=0, vary=False)  # Silicate ratio
ptype1.add('pt1_U_Si', value=2.0, vary=False)   # Silicate Half Saturation Constant
ptype1.add('pt1_U_N', value=0.446, vary=False)    # Nitrate Half Saturation Constant
ptype1.add('pt1_muP', value=1.5, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 2 (e.g. Coccos)
ptype2 = Parameters()
ptype2.add('pt2_U_N', value=0.265, vary=False)    # Nitrate Half Saturation Constant
ptype2.add('pt2_muP', value=1.1, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 3 (e.g. Dinos)
ptype3 = Parameters()
ptype3.add('pt3_U_N', value=0.009, vary=False)    # Nitrate Half Saturation Constant
ptype3.add('pt3_muP', value=0.5, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 3 (e.g. Nanos)
ptype4 = Parameters()
ptype4.add('pt4_U_N', value=0.045, vary=False)    # Nitrate Half Saturation Constant
ptype4.add('pt4_muP', value=1.7, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# z - related
#z grazing related
standardparams.add('moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
standardparams.add('deltaZ', value=0.75, vary=False)    # Zooplankton Grazing assimilation coefficient (-)
standardparams.add('deltaLambda', value=0.75, vary=False)    # Zooplankton Inter-Grazing assimilation coefficient (-)
standardparams.add('muIntGraze', value=0.01, vary=False)  # InterZooGrazing maximum grazing rate
standardparams.add('kIntGraze', value=0.5, vary=False)  # InterZooGrazing saturation constant

standardparams.add('Kp', value=0, vary=False)     # Zooplankton Grazing saturation constant (-)
standardparams.add('pred', value=0, vary=False)  # quadratic higher order predation rate on zooplankton
standardparams.add('muZ', value=0, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# set up zooplankton type 1 (e.g. MIKRO zooplankton)
ztype1 = Parameters()
ztype1.add('zt1_muZ', value=0.1, vary=False)    # Zooplankton maximum grazing rate (d^-1)

ztype1.add('zt1_Kp', value=0.5, vary=False)       # Zooplankton Grazing saturation constant (-)
ztype1.add('zt1_pred', value=0.01, vary=False)    # quadratic higher order predation rate on zooplankton

# set up zooplankton type 2 (e.g. MESO zooplankton)
ztype2 = Parameters()
ztype2.add('zt2_muZ', value=0.1, vary=False)    # Zooplankton maximum grazing rate (d^-1)

ztype2.add('zt2_Kp', value=0.5, vary=False)       # Zooplankton Grazing saturation constant (-)
ztype2.add('zt2_pred', value=0.01, vary=False)    # quadratic higher order predation rate on zooplankton

"""
add feeding pref params!
Z1P1
Z1P2
Z2P1
Z2P2
"""
# MIKRO
ztype1.add('zt1_P1', value=0, vary=False)  # Diatoms
ztype1.add('zt1_P2', value=1, vary=False)  # Coccos
ztype1.add('zt1_P3', value=1, vary=False)  # Dinos
ztype1.add('zt1_P4', value=1, vary=False)  # Nano
# MESO
ztype2.add('zt2_P1', value=1, vary=False)
ztype2.add('zt2_P2', value=1, vary=False)
ztype2.add('zt2_P3', value=1, vary=False)
ztype2.add('zt2_P4', value=0, vary=False)

ptype1.add('pt1_Z1', value=ztype1['zt1_P1'].value, vary=False)
ptype1.add('pt1_Z2', value=ztype2['zt2_P1'].value, vary=False)

ptype2.add('pt2_Z1', value=ztype1['zt1_P2'].value, vary=False)
ptype2.add('pt2_Z2', value=ztype2['zt2_P2'].value, vary=False)

ptype3.add('pt3_Z1', value=ztype1['zt1_P3'].value, vary=False)
ptype3.add('pt3_Z2', value=ztype2['zt2_P3'].value, vary=False)

ptype4.add('pt4_Z1', value=ztype1['zt1_P4'].value, vary=False)
ptype4.add('pt4_Z2', value=ztype2['zt2_P4'].value, vary=False)


def setupinitcond(pfn,zn):
    # initialize parameters:
    N0 = np.mean(mc.NOX)  # Initial Nitrate concentration (mmol*m^-3)
    Si0 = np.mean(mc.SiOX)  # Initial Silicate concentration (mmol*m^-3)
    Z0 = 0.1 / zn  # Initial Zooplankton concentration (mmol*m^-3)
    D0 = 0.01  # Initial Detritus concentration (mmol*m^-3)
    P0 = 0.01 / pfn  # Initial Phytoplankton concentration (mmol*m^-3)

    initnut = [N0, Si0, D0]
    initzoo = [Z0 for i in range(zn)]
    initphy = [P0 for i in range(pfn)]
    outputl = [0 for i in range(20)]
    initcond = np.concatenate([initnut, initzoo, initphy, outputl])
    #print(type(initcond))
    return initcond


def runmodel(all_params, initcond):
    print(list(all_params)[:])

    z = Plankton(all_params, 'Zooplankton').init()
    p = Plankton(all_params, 'Phytoplankton').init()

    # INTEGRATE:
    tos = time.time()
    print('starting integration')
    outarray = odeint(mc.phytomftm_extendedoutput, initcond, timedays_model, args=(all_params, p, z))#, rtol=1e-12, atol=1e-12)
    tos1 = time.time()
    print('finished after %4.3f sec' % (tos1 - tos))

    return outarray


def callmodelrun(pfn,zn):
    # number of phytoplankton func types
    standardparams.add('pfun_num', value=pfn, vary=False)
    # number of zooplankton groups
    standardparams.add('zoo_num', value=zn, vary=False)

    if pfn == 4 and zn == 2:
         print('4P2Z - prelim model')
         all_params = (standardparams + ptype1 + ptype2 + ptype3 + ptype4 + ztype1 + ztype2)
    elif pfn == 2 and zn == 2:
         print('2P2Z')
         all_params = (standardparams + ptype1 + ptype2 + ztype1 + ztype2)
    elif pfn == 2 and zn == 1:
         print('2P1Z')
         all_params = (standardparams + ptype1 + ptype2 + ztype1)
    elif pfn == 1 and zn == 2:
         print('2P1Z')
         all_params = (standardparams + ptype1 + ztype1 + ztype2)
    elif pfn == 1 and zn == 1:
         print('1P1Z')
         all_params = (standardparams + ptype1 + ztype1)
    else:
        print('just standard params')
        all_params = (standardparams)

    #all_params = (standardparams)
    parameters = all_params
    initialcond = setupinitcond(pfn,zn)
    print(initialcond)
    out = runmodel(parameters,initialcond)

    return out


def plotoutput(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays = timedays_model#[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray#[1460:1825]

    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    ax1[i_plot].set_title(title)

    # Figure 1
    # N
    ax1[i_plot].plot(timedays, outarray_ly[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
    if i_plot == 0:
        ax1[i_plot].set_ylabel('Nitrate \n' '[µM]', multialignment='center', fontsize=10)
    #ax1[i_plot].set_ylim(-0.1, 5)

    # Si
    ax2[i_plot].plot(timedays, outarray_ly[:, 1], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax2[i_plot].set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    #ax2[i_plot].set_ylim(-0.1, 12)

    # Phyto
    ax3[i_plot].plot(timedays, sum([outarray_ly[:, 3 + zn + i] for i in range(pfn)]), c=colors[4], lw=lws[1])
    [ax3[i_plot].plot(timedays, outarray_ly[:, 3 + zn + i], c=colors[i + 1]) for i in range(pfn)]
    if i_plot == 0:
        ax3[i_plot].set_ylabel('Phyto \n' '[µM N]', multialignment='center', fontsize=10)
    #ax3[i_plot].set_ylim(-0.1, 0.8)

    #ax3[i_plot].set_title('Phy Biomass & ChlA Data')

    # Z
    ax4[i_plot].plot(timedays, sum([outarray_ly[:, 3 + i] for i in range(zn)]), c=colors[4], lw=lws[1])
    [ax4[i_plot].plot(timedays, outarray_ly[:, 3 + i], c=colors[i + 1], lw=lws[0], alpha=alphas[0]) for i in range(zn)]
    if i_plot == 0:
        ax4[i_plot].set_ylabel('Zooplankton \n' '[µM N]', multialignment='center', fontsize=9)
    ax4[i_plot].tick_params('y', labelsize=10)

    #ax4[i_plot].set_title('Zooplankton')
    #ax4[i_plot].set_ylim(0, 0.62)

    # D
    ax5[i_plot].plot(timedays, outarray_ly[:, 2], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax5[i_plot].set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    #ax5[i_plot].set_title('Detritus')
    #ax5[i_plot].set_ylim(0,0.15)

    ax5[i_plot].set_xlabel('Day in year')
    # Legend

def plot1output(outarray, pfn, zn, i_plot, title):
    f1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex='col', sharey='row')

    # PLOTTING
    timedays = timedays_model#[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray#[1460:1825]

    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    ax1.set_title(title)

    # Figure 1
    # N
    ax1.plot(timedays, outarray_ly[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
    if i_plot == 0:
        ax1.set_ylabel('Nitrate \n' '[µM]', multialignment='center', fontsize=10)
    #ax1[i_plot].set_ylim(-0.1, 5)

    # Si
    ax2.plot(timedays, outarray_ly[:, 1], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax2.set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    #ax2[i_plot].set_ylim(-0.1, 12)

    # Phyto
    ax3.plot(timedays, sum([outarray_ly[:, 3 + zn + i] for i in range(pfn)]), c=colors[4], lw=lws[1])
    [ax3.plot(timedays, outarray_ly[:, 3 + zn + i], c=colors[i + 1]) for i in range(pfn)]
    if i_plot == 0:
        ax3.set_ylabel('Phyto \n' '[µM N]', multialignment='center', fontsize=10)
    #ax3[i_plot].set_ylim(-0.1, 0.8)

    #ax3[i_plot].set_title('Phy Biomass & ChlA Data')

    # Z
    ax4.plot(timedays, sum([outarray_ly[:, 3 + i] for i in range(zn)]), c=colors[4], lw=lws[1])
    [ax4.plot(timedays, outarray_ly[:, 3 + i], c=colors[i + 1], lw=lws[0], alpha=alphas[0]) for i in range(zn)]
    if i_plot == 0:
        ax4.set_ylabel('Zooplankton \n' '[µM N]', multialignment='center', fontsize=9)
    ax4.tick_params('y', labelsize=10)

    #ax4[i_plot].set_title('Zooplankton')
    #ax4[i_plot].set_ylim(0, 0.62)

    # D
    ax5.plot(timedays, outarray_ly[:, 2], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax5.set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    #ax5[i_plot].set_title('Detritus')
    #ax5[i_plot].set_ylim(0,0.15)

    ax5.set_xlabel('Day in year')
    # Legend

def plotYEARoutput(outarray, pfn, zn, i_plot, title):
    f1, (ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9) = plt.subplots(9, 1, sharex='col', sharey='row')
    # N / Si / Pdt / Pc / Pdn / Pn / Zmu / Zlambda / D
    # PLOTTING
    timedays = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]
    i_plot=0
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [1, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

    ax1.set_title(title)

    # Figure 1
    # N
    ax1.plot(timedays, outarray_ly[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
    if i_plot == 0:
        ax1.set_ylabel('Nitrate \n' '[µM]', multialignment='center', fontsize=10)
    #ax1[i_plot].set_ylim(-0.1, 5)

    # Si
    ax2.plot(timedays, outarray_ly[:, 1], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax2.set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    #ax2[i_plot].set_ylim(-0.1, 12)

    # Phyto
    i=2
    #ax3.plot(timedays, sum([outarray_ly[:, 3 + zn + i] for i in range(pfn)]), c=colors[4], lw=lws[1])
    ax3.plot(timedays, outarray_ly[:, 3 + zn + 0], c=colors[i + 1])
    ax4.plot(timedays, outarray_ly[:, 3 + zn + 1], c=colors[i + 1])
    ax5.plot(timedays, outarray_ly[:, 3 + zn + 2], c=colors[i + 1])
    ax6.plot(timedays, outarray_ly[:, 3 + zn + 3], c=colors[i + 1])

    if i_plot == 0:
        ax3.set_ylabel('Diatom \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax4.set_ylabel('Coccs \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax5.set_ylabel('Dinos \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax6.set_ylabel('Nano \n' '[µM N]', multialignment='center', fontsize=10)
    ax3.set_ylim(-0.1, 0.8)
    ax4.set_ylim(-0.1, 0.8)
    ax5.set_ylim(-0.1, 0.8)
    ax6.set_ylim(-0.1, 0.8)

    #ax3[i_plot].set_title('Phy Biomass & ChlA Data')

    # Z
    i=3
    #ax4.plot(timedays, sum([outarray_ly[:, 3 + i] for i in range(zn)]), c=colors[4], lw=lws[1])
    #[ax4.plot(timedays, outarray_ly[:, 3 + i], c=colors[i + 1], lw=lws[0], alpha=alphas[0]) for i in range(zn)]
    ax7.plot(timedays, outarray_ly[:, 3 + 0], c=colors[i + 1], lw=lws[0], alpha=alphas[0])
    ax8.plot(timedays, outarray_ly[:, 3 + 1], c=colors[i + 1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax7.set_ylabel('Mikro Z \n' '[µM N]', multialignment='center', fontsize=9)
    if i_plot == 0:
        ax8.set_ylabel('Meso Z \n' '[µM N]', multialignment='center', fontsize=9)
    ax7.tick_params('y', labelsize=10)

    #ax4[i_plot].set_title('Zooplankton')
    #ax4[i_plot].set_ylim(0, 0.62)

    # D
    ax9.plot(timedays, outarray_ly[:, 2], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax9.set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    #ax5[i_plot].set_title('Detritus')
    #ax5[i_plot].set_ylim(0,0.15)

    ax9.set_xlabel('Day in year')
    # Legend


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
#FFDE00 Silicate
#575756 Nitrate
#009640 Phytoplankton
#EA5B0C Mikrozooplankton
#BE1622 Mesozooplankton
#B17F4A Detritus
#DEDC00 SST
#F9B233 PAR

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
    i_plot=0
    # color vectors
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
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
    #ax1[i_plot].set_ylim(-0.1, 5)

    # Si
    ax2.plot(timedays, outarray_ly[:, 1], c="Si", lw=lws[0], alpha=alphas[0])
    # Si Data
    ax2.scatter(SiOH['yday'].values, SiOH['SiOH'].values, c="Si", s=4.3)
    if i_plot == 0:
        ax2.set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    #ax2[i_plot].set_ylim(-0.1, 12)

    # Phyto
    i=2
    #ax3.plot(timedays, sum([outarray_ly[:, 3 + zn + i] for i in range(pfn)]), c=colors[4], lw=lws[1])
    ax3.plot(timedays, outarray_ly[:, 3 + zn + 0], c=colors[i + 1])
    ax4.plot(timedays, outarray_ly[:, 3 + zn + 1], c=colors[i + 1])
    ax5.plot(timedays, outarray_ly[:, 3 + zn + 2], c=colors[i + 1])
    ax6.plot(timedays, outarray_ly[:, 3 + zn + 3], c=colors[i + 1])

    if i_plot == 0:
        ax3.set_ylabel('Diatom \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax4.set_ylabel('Coccs \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax5.set_ylabel('Dinos \n' '[µM N]', multialignment='center', fontsize=10)
    if i_plot == 0:
        ax6.set_ylabel('Nano \n' '[µM N]', multialignment='center', fontsize=10)
    ax3.set_ylim(-0.1, 0.8)
    ax4.set_ylim(-0.1, 0.8)
    ax5.set_ylim(-0.1, 0.8)
    ax6.set_ylim(-0.1, 0.8)

    #ax3[i_plot].set_title('Phy Biomass & ChlA Data')

    # Z
    i=3
    #ax4.plot(timedays, sum([outarray_ly[:, 3 + i] for i in range(zn)]), c=colors[4], lw=lws[1])
    #[ax4.plot(timedays, outarray_ly[:, 3 + i], c=colors[i + 1], lw=lws[0], alpha=alphas[0]) for i in range(zn)]
    ax7.plot(timedays, outarray_ly[:, 3 + 0], c=colors[i + 1], lw=lws[0], alpha=alphas[0])
    ax8.plot(timedays, outarray_ly[:, 3 + 1], c=colors[i + 1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax7.set_ylabel('Mikro Z \n' '[µM N]', multialignment='center', fontsize=9)
    if i_plot == 0:
        ax8.set_ylabel('Meso Z \n' '[µM N]', multialignment='center', fontsize=9)
    ax7.tick_params('y', labelsize=10)

    #ax4[i_plot].set_title('Zooplankton')
    #ax4[i_plot].set_ylim(0, 0.62)

    # D
    ax9.plot(timedays, outarray_ly[:, 2], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax9.set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    #ax5[i_plot].set_title('Detritus')
    #ax5[i_plot].set_ylim(0,0.15)

    ax9.set_xlabel('Day in year')
    # Legend

timedays_model = np.arange(0., 5 * 365., 1.0)

out4P2Z = callmodelrun(4,2)

plotDATAvsYEARoutput(out4P2Z, 4, 2, 1, '4P2Z')

#out2P2Z = callmodelrun(2,2)

#plot1output(out2P2Z, 2, 2, 1, '1P1Z')