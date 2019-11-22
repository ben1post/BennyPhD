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
from PhytoMFTM.ModelClasses import Plankton, Forcing
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
standardparams.add('kappa', value=0.1, vary=False) # vary=False) # min=0.09, max=0.11) # Diffusive mixing across thermocline (m*d^-1)
standardparams.add('deltaD_N', value=0.1, vary=False)   # Nitrate Remineralization rate (d^-1)

standardparams.add('kw', value=0.04, vary=False)     # Light attenuation constant of water (m^-1)

standardparams.add('kc', value=0.03, vary=False)      # Light attenuation via phytoplankton pigment (m^-1)
standardparams.add('alpha', value=0.15, vary=False)  # initial slope of the P-I curve
standardparams.add('VpMax', value=1., vary=False)    # maximum photosynthetic rate

standardparams.add('v', value=0.01, vary=False)      # Sinking of Phytoplankton from Mixed Layer
standardparams.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)

# p - related
standardparams.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)

standardparams.add('ratioSi', value=0, vary=False)  # Silicate ratio ## FALLBACK PARAM FOR OTHER PFTs
standardparams.add('U_N', value=0, vary=False)    # Nitrate Half Saturation Constant
standardparams.add('U_Si', value=0, vary=False)   # Silicate Half Saturation Constant
standardparams.add('muP', value=0, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# set up phytoplankton type 1 (e.g. DIATOMS)
ptype1 = Parameters()
ptype1.add('pt1_ratioSi', value=0., vary=False)  # Silicate ratio
ptype1.add('pt1_U_Si', value=0., vary=False)   # Silicate Half Saturation Constant
ptype1.add('pt1_U_N', value=0.9, vary=False)    # Nitrate Half Saturation Constant
ptype1.add('pt1_muP', value=1., vary=False)    # Phytoplankton maximum growth rate (d^-1)


# z - related
#z grazing related
standardparams.add('moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
standardparams.add('deltaZ', value=0.75, vary=False)    # Zooplankton Grazing assimilation coefficient (-)
standardparams.add('deltaLambda', value=0.75, vary=False)    # Zooplankton Inter-Grazing assimilation coefficient (-)
standardparams.add('muIntGraze', value=0., vary=False)  # InterZooGrazing maximum grazing rate
standardparams.add('kIntGraze', value=0., vary=False)  # InterZooGrazing saturation constant

standardparams.add('Kp', value=0, vary=False)     # Zooplankton Grazing saturation constant (-)
standardparams.add('pred', value=0, vary=False)  # quadratic higher order predation rate on zooplankton
standardparams.add('muZ', value=0, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# set up zooplankton type 1 (e.g. MIKRO zooplankton)
ztype1 = Parameters()
ztype1.add('zt1_muZ', value=1.,  vary=False)    # Zooplankton maximum grazing rate (d^-1)

ztype1.add('zt1_Kp', value=2.5,  vary=False)       # Zooplankton Grazing saturation constant (-)
ztype1.add('zt1_pred', value=0.01,  vary=False)    # quadratic higher order predation rate on zooplankton


# MIKRO
ztype1.add('zt1_P1', value=1, vary=False)  # Diatoms

# MESO

# inter zoo feeding
# MIKRO
ztype1.add('zt1_Zint_feed1', value=0, vary=False)
ztype1.add('zt1_Zint_feed2', value=0, vary=False)

# CONVERT FEEDPREFS TO GRAZEPREF FOR CALCULATION OF GRAZING
ztype1.add('zt1_Zint_grazed1', value=ztype1['zt1_Zint_feed1'].value, vary=False)

ptype1.add('pt1_Z1', value=ztype1['zt1_P1'].value, vary=False)


#BIOTRANS_ChlA = pandas.read_csv('ValidationData/BIOTRANS_SeaWIFs_8Day_ChlA_extracted.csv')

#PAPA_ChlA = pandas.read_csv('ValidationData/PAPA_SeaWIFs_8Day_ChlA_extracted.csv')


# set up model conditions and parameter dict
def setupinitcond(pfn,zn):
    # initialize parameters:
    N0 = 2  # Initial Nitrate concentration (mmol*m^-3)
    Si0 = 0  # Initial Silicate concentration (mmol*m^-3)
    Z0 = 0.1   # Initial Zooplankton concentration (mmol*m^-3)
    D0 = 0.1  # Initial Detritus concentration (mmol*m^-3)
    P0 = 1.0  # Initial Phytoplankton concentration (mmol*m^-3)

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
    outarray = odeint(mc.phytomftm_extendedoutput_forcing, initcond, timedays_model, args=(all_params, p, z, fx))
    tos1 = time.time()
    print('finished after %4.3f sec' % (tos1 - tos))

    return outarray


def callmodelrun(pfn,zn):
    # number of phytoplankton func types
    standardparams.add('pfun_num', value=pfn, vary=False)
    # number of zooplankton groups
    standardparams.add('zoo_num', value=zn, vary=False)
    """
    if pfn == 2 and zn == 2:
         print('3P2Z')
         all_params = (standardparams + ztype1 + ztype2)
    elif pfn == 3 and zn == 2:
         print('2P2Z')
         all_params = (standardparams + ztype1 + ztype2)
    elif pfn == 4 and zn == 1:
         print('2P1Z')
         all_params = (standardparams )
    elif pfn == 1 and zn == 1:
         print('1P1Z')
         all_params = (standardparams)
    else:
        print('just standard params')
        all_params = (standardparams)
    """
    all_params = (standardparams + ptype1 + ztype1)

    parameters = all_params
    initialcond = setupinitcond(pfn,zn)
    print(initialcond)
    out = runmodel(parameters,initialcond)

    return out

def plotoutput(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays = timedays_model[1:366]#[1:120]#[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]#[1:120]#[1460:1825]

    # color vectors
    #colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', 'grey']
    alphas = [1., 0.8, 0.6, 0.4]
    lws = [2, 2.5, 4, 5.5]

    # artist for legends
    FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])


    ax1[i_plot].set_title(title)

    # Figure 1
    # N
    N_Max = np.max(fx.NOX.return_interpvalattime(timedays)) + np.max(fx.NOX.return_interpvalattime(timedays)) * 0.1
    ax1[i_plot].plot(timedays, outarray_ly[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
    if i_plot == 0:
        ax1[i_plot].set_ylabel('Nutrients \n' '[µM N]', multialignment='center', fontsize=10)
    ax1[i_plot].set_ylim(0, N_Max)

    # Si
    #ax2[i_plot].plot(timedays, outarray_ly[:, 1], c=colors[1], lw=lws[0], alpha=alphas[0])
    #if i_plot == 0:
    #    ax2[i_plot].set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    #ax2[i_plot].set_ylim(-0.1, 12)

    Chl_Max = 1
    # Phyto
    ax2s = ax2[i_plot].twinx()
    ax2s.scatter(np.arange(1,13,1) * 30 -15, fx.chla.MODISaquaChlAclimatology)
    ax2[i_plot].set_zorder(ax2s.get_zorder() + 1)
    ax2[i_plot].patch.set_visible(False)
    ax2s.set_ylim(0, Chl_Max)
    ax2s.set_ylabel('$Chl$ $a$ \n' '[mg $m^{-3}$]', multialignment='center', fontsize=10)
    # UNITS: mg m^-3 according to ncdf

    Pall = [outarray_ly[:, 3 + zn + i] for i in range(pfn)]
    P_Max = np.max(Pall) + 0.5 * np.max(Pall)


    ax2[i_plot].plot(timedays, sum([outarray_ly[:, 3 + zn + i] for i in range(pfn)]), c=colors[4], lw=lws[1])
    [ax2[i_plot].plot(timedays, outarray_ly[:, 3 + zn + i], c=colors[i + 1]) for i in range(pfn)]
    if i_plot == 0:
        ax2[i_plot].set_ylabel('Phytoplankton \n' '[µM N]', multialignment='center', fontsize=10)
    ax2[i_plot].set_ylim(0, P_Max)

    #ax3[i_plot].set_title('Phy Biomass & ChlA Data')

    # Z
    Zall = [outarray_ly[:, 3 + i] for i in range(zn)]
    Z_Max = np.max(Zall) + 0.1 * np.max(Zall)

    ax3[i_plot].plot(timedays, sum([outarray_ly[:, 3 + i] for i in range(zn)]), c=colors[4], lw=lws[1])
    [ax3[i_plot].plot(timedays, outarray_ly[:, 3 + i], c=colors[i + 1], lw=lws[0], alpha=alphas[0]) for i in range(zn)]
    if i_plot == 0:
        ax3[i_plot].set_ylabel('Zooplankton \n' '[µM N]', multialignment='center', fontsize=9)
    ax3[i_plot].tick_params('y', labelsize=10)
    ax3[i_plot].set_ylim(0, Z_Max)
    #ax4[i_plot].set_title('Zooplankton')

    D_Max = np.max(outarray_ly[:, 2]) + 0.2 * np.max(outarray_ly[:, 2])
    # D
    ax4[i_plot].plot(timedays, outarray_ly[:, 2], c=colors[1], lw=lws[0], alpha=alphas[0])
    if i_plot == 0:
        ax4[i_plot].set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    #ax5[i_plot].set_title('Detritus')
    ax4[i_plot].set_ylim(0,D_Max)

    ax4[i_plot].set_xlabel('Day in year')
    # Legend

    ## PHYSICS ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    muplot = i_plot+1
    ax1[muplot].set_title('model forcing')

    ax4[muplot].set_xlabel('Day in year')

    ax1[muplot].plot(timedays, fx.NOX.return_interpvalattime(timedays), c=colors[5], lw=lws[0], alpha=alphas[0])
    ax1[muplot].set_ylabel('$N_0$ \n' '[µmol $kg^{-1}$]', multialignment='center', fontsize=10)
    ax1[muplot].set_ylim(0., N_Max)
    #ax1[muplot].invert_yaxis()

    MLD_max = np.max(fx.MLD.return_interpvalattime(timedays)) + 0.1 * np.max(fx.MLD.return_interpvalattime(timedays))
    ax2[muplot].plot(timedays, fx.MLD.return_interpvalattime(timedays), c=colors[5], lw=lws[0], alpha=alphas[0])
    ax2[muplot].set_ylabel('MLD \n' '[m]', multialignment='center', fontsize=10)
    ax2[muplot].set_ylim(0, MLD_max) # 400 for biotrans, 100 for Papa
    ax2[muplot].invert_yaxis()

    PAR_max = np.max(fx.PAR.return_interpvalattime(timedays)) + 0.1 * np.max(fx.PAR.return_interpvalattime(timedays))
    ax3[muplot].plot(timedays, fx.PAR.return_interpvalattime(timedays), c=colors[5], lw=lws[0], alpha=alphas[0])
    ax3[muplot].set_ylabel('PAR \n' '[E $m^{−2}$ $s^{−1}$]', multialignment='center', fontsize=10)
    ax3[muplot].set_ylim(0, PAR_max)
    # ax1[muplot].invert_yaxis()

    Tmld_max = np.max(fx.SST.return_interpvalattime(timedays)) + 0.1 * np.max(fx.SST.return_interpvalattime(timedays))
    ax4[muplot].plot(timedays, fx.SST.return_interpvalattime(timedays), c=colors[5], lw=lws[0], alpha=alphas[0])
    ax4[muplot].set_ylabel('$T_{MLD}$ \n' '[°C]', multialignment='center', fontsize=10)
    ax4[muplot].set_ylim(0, Tmld_max)
    # ax1[muplot].invert_yaxis()








timedays_model = np.arange(0., 5 * 365., 1.0)

numcols = 2

fx = Forcing('PAPA')

out1P1Z = callmodelrun(1,1)


f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, numcols, sharex='col')#, sharey='row')


plt.setp((ax1, ax2, ax3, ax4), xticks=[1,60,120,180,240,300,365])
from matplotlib.ticker import MaxNLocator
for axe in (ax1, ax2, ax3, ax4):
    for i in range(numcols):
        axe[i].get_yaxis().set_major_locator(MaxNLocator(nbins=4))

# remove labels from axes:
#for axe in (ax1, ax2, ax3, ax4):
#    for i in [0, 1]:
#        axe[i].yaxis.set_major_formatter(plt.NullFormatter())



plotoutput(out1P1Z,1,1,0, '1P1Z')

plt.tight_layout()