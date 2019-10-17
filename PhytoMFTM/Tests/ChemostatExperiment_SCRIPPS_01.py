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

# set up basic parameters of model:
standardparams = Parameters()

# NEW PARAM: 'flow',

# UNNECESSARY PARAMS:
# mld - related
standardparams.add('kappa', value=0, vary=False) # vary=False) # min=0.09, max=0.11) # Diffusive mixing across thermocline (m*d^-1)
standardparams.add('kw', value=0.04, vary=False)     # Light attenuation constant of water (m^-1)

standardparams.add('kc', value=0.03, vary=False)      # Light attenuation via phytoplankton pigment (m^-1)
standardparams.add('alpha', value=0.15, vary=False)  # initial slope of the P-I curve
standardparams.add('VpMax', value=1., vary=False)    # maximum photosynthetic rate

standardparams.add('v', value=0., vary=False)      # Sinking of Phytoplankton from Mixed Layer
standardparams.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)



# NECESSARY PARAMS:
standardparams.add('flow', value=0.1, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)
standardparams.add('deltaD_N', value=0, vary=False)   # Nitrate Remineralization rate (d^-1)

# p - related
standardparams.add('moP', value=0.01, vary=False)    # Phytoplankton mortality (d^-1)

standardparams.add('ratioSi', value=0, vary=False)  # Silicate ratio ## FALLBACK PARAM FOR OTHER PFTs
standardparams.add('U_N', value=0, vary=False)    # Nitrate Half Saturation Constant
standardparams.add('U_Si', value=0, vary=False)   # Silicate Half Saturation Constant
standardparams.add('muP', value=0, vary=False)    # Phytoplankton maximum growth rate (d^-1)

# z - related
#z grazing related
standardparams.add('moZ', value=0.07, vary=False)        # Zooplankton mortality (d^-1)
standardparams.add('deltaZ', value=0.75, vary=False)    # Zooplankton Grazing assimilation coefficient (-)
standardparams.add('deltaLambda', value=0.75, vary=False)    # Zooplankton Inter-Grazing assimilation coefficient (-)
standardparams.add('muIntGraze', value=0.1, vary=False)  # InterZooGrazing maximum grazing rate
standardparams.add('kIntGraze', value=0.5, vary=False)  # InterZooGrazing saturation constant

standardparams.add('Kp', value=0, vary=False)     # Zooplankton Grazing saturation constant (-)
standardparams.add('pred', value=0, vary=False)  # quadratic higher order predation rate on zooplankton
standardparams.add('muZ', value=0, vary=False)    # Zooplankton maximum grazing rate (d^-1)


# set up zooplankton type 1 (e.g. MIKRO zooplankton)
ztype1 = Parameters()
ztype1.add('zt1_muZ', value=0.5, min=.2, max=1.5)    # Zooplankton maximum grazing rate (d^-1)

ztype1.add('zt1_Kp', value=1., min=.2, max=1.5)       # Zooplankton Grazing saturation constant (-)
ztype1.add('zt1_pred', value=0.06, min=.001, max=0.2)    # quadratic higher order predation rate on zooplankton


"""
add feeding pref params!
Z1P1
Z1P2
Z2P1
Z2P2
"""
'''''PARAM GENERATOR'''
# here!



# set up model conditions and parameter dict
def setupinitcond(pfn,zn):
    # initialize parameters:
    N0 = 1  # Initial Nitrate concentration (mmol*m^-3)
    Si0 = 0  # Initial Silicate concentration (mmol*m^-3)
    Z0 = 0.   # Initial Zooplankton concentration (mmol*m^-3)
    D0 = 0.0  # Initial Detritus concentration (mmol*m^-3)
    P0 = 0.5 / pfn  # Initial Phytoplankton concentration (mmol*m^-3)

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
    outarray = odeint(mc.phytomftm_chemostat, initcond, timedays_model, args=(all_params, p, z, fx))
    tos1 = time.time()
    print('finished after %4.3f sec' % (tos1 - tos))

    return outarray


def callmodelrun(pfn,zn):
    # number of phytoplankton func types
    standardparams.add('pfun_num', value=pfn, vary=False)
    # number of zooplankton groups
    standardparams.add('zoo_num', value=zn, vary=False)

    standardparams.add('Zint_feed1', value=0, vary=False)  # Phytoplankton maximum growth rate (d^-1)
    standardparams.add('Zint_grazed1', value=0, vary=False)
    standardparams.add('Z1', value=1, vary=False)

    # uptake * growth = constant
    # 1. give mean value for both
    # 2. generate ranges
    mean_U_N = 1
    range_U_N = 0.5

    if pfn == 1:
        uptakeRange = mean_U_N

        growthRange = 1 / uptakeRange

        ptype = Parameters()

        ptype.add('pt1_U_N', value=uptakeRange, vary=False)
        ptype.add('pt1_muP', value=growthRange, vary=False)
        ptype.add('P1', value=1, vary=False)

        all_params = standardparams + ztype1 + ptype

    else:
        uptakeRange = np.linspace(mean_U_N - range_U_N, mean_U_N + range_U_N, pfn)

        growthRange = 1 / uptakeRange

        ptype = Parameters()
        for n in range(pfn):
            ptype.add('pt' + str(n + 1) + '_U_N', value=uptakeRange[n], vary=False)  # uptake affinity
            ptype.add('pt' + str(n + 1) + '_muP', value=growthRange[n], vary=False)  # growth rate
            ptype.add('P' + str(n + 1), value=1 / pfn, vary=False)  # edibility

        print(ptype)
        all_params = standardparams + ztype1 + ptype





    parameters = all_params
    initialcond = setupinitcond(pfn,zn)
    print(initialcond)
    out = runmodel(parameters,initialcond)

    return out

def plotoutput(outarray, pfn, zn, i_plot, title):

    # PLOTTING
    timedays = timedays_model[1:120]#[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1:120]#[1460:1825]

    # color vectors
    #colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
    colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73','#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73','#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', '#009e73']
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
    #ax2[i_plot].plot(timedays, outarray_ly[:, 1], c=colors[1], lw=lws[0], alpha=alphas[0])
    #if i_plot == 0:
    #    ax2[i_plot].set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
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


""" here the functions end, and the code calling them starts!! """



timedays_model = np.arange(0., 5 * 365., 1.0)
fx = Forcing('flowthrough')


out1P1Z = callmodelrun(1, 1)
out5P1Z = callmodelrun(5, 1)
out15P1Z = callmodelrun(15, 1)


f1, (ax1, ax3, ax4, ax5) = plt.subplots(4, 3, sharex='col', sharey='row') # ax2 removed for single nut plots
#plt.subplots_adjust(hspace=0.01)

plotoutput(out1P1Z,1,1,0, '1P1Z')

plotoutput(out5P1Z,5,1,1, '5P1Z')

plotoutput(out15P1Z,15,1,2, '15P1Z')

