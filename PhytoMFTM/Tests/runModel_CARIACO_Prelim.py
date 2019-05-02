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
standardparams.add('deltaD_Si', value=0.05, vary=False)  # Silicate Mineralization rate (d^-1)

# z - related
standardparams.add('moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
standardparams.add('deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)
#z grazing related
standardparams.add('gr_p', value=0.6, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
standardparams.add('muZ', value=0.4, vary=False)    # Zooplankton maximum grazing rate (d^-1)

# p - related
standardparams.add('kw', value=0.1, vary=False)     # Light attenuation constant (m^-1)
standardparams.add('OptI', value=30, vary=False)    # Optimum irradiance (einstein*m^-2*d^-1)

standardparams.add('U_N', value=5.0, vary=False)    # Nitrate Half Saturation Constant
standardparams.add('U_Si', value=5.0, vary=False)   # Silicate Half Saturation Constant

standardparams.add('v', value=0.11, vary=False)      # Sinking of Phytoplankton from Mixed Layer

standardparams.add('muP', value=1.6, vary=False)    # Phytoplankton maximum growth rate (d^-1)
standardparams.add('moP', value=0.1, vary=False)    # Phytoplankton mortality (d^-1)

standardparams.add('Kp', value=0.3, vary=False)     # Zooplankton Grazing saturation constant (-)
standardparams.add('ratioSi', value=1.1, vary=False)  # Silicate ratio

# set up phytoplankton type 1 (e.g. Nanoflagellates)
ptype1 = Parameters()
ptype1.add('pt1_U_Si', value=5.0, vary=False)   # Silicate Half Saturation Constant
ptype1.add('pt1_muP', value=1.4, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype1.add('pt1_ratioSi', value=1.25, vary=False)  # Silicate ratio

# set up phytoplankton type 2 (e.g. Diatoms)
ptype2 = Parameters()
ptype2.add('pt2_U_Si', value=0., vary=False)     # Silicate Half Saturation Constant
ptype2.add('pt2_muP', value=1.1, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype2.add('pt2_ratioSi', value=0., vary=False)  # Silicate ratio

# set up phytoplankton type 3 (e.g. Coccos)
ptype3 = Parameters()
ptype3.add('pt3_U_Si', value=0., vary=False)     # Silicate Half Saturation Constant
ptype3.add('pt3_U_N', value=1.94, vary=False)    # Nitrate Half Saturation Constant
ptype3.add('pt3_muP', value=0.5, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype3.add('pt3_ratioSi', value=0., vary=False)  # Silicate ratio

# set up phytoplankton type 3 (e.g. Dinos)
ptype4 = Parameters()
ptype4.add('pt3_U_Si', value=0., vary=False)     # Silicate Half Saturation Constant
ptype4.add('pt3_U_N', value=1.94, vary=False)    # Nitrate Half Saturation Constant
ptype4.add('pt3_muP', value=0.5, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype4.add('pt3_ratioSi', value=0., vary=False)  # Silicate ratio


# set up zooplankton type 1 (e.g. small zooplankton)
ztype1 = Parameters()
ztype1.add('zt1_gr_p', value=0.6, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
ztype1.add('zt1_muZ', value=0.05, vary=False)    # Zooplankton maximum grazing rate (d^-1)

ztype1.add('zt1_moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
ztype1.add('zt1_deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)

# set up zooplankton type 2 (e.g. larger zooplankton)
ztype2 = Parameters()
ztype2.add('zt2_gr_p', value=0.6, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
ztype2.add('zt2_muZ', value=0.05, vary=False)    # Zooplankton maximum grazing rate (d^-1)

ztype2.add('zt2_moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
ztype2.add('zt2_deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)

"""
add feeding pref params!
Z1P1
Z1P2
Z2P1
Z2P2
"""
ztype1.add('zt1_P1', value=1, vary=False)
ztype1.add('zt1_P2', value=.2, vary=False)
ztype1.add('zt1_P3', value=.3, vary=False)
ztype1.add('zt1_P4', value=.4, vary=False)

ztype2.add('zt2_P1', value=1.1, vary=False)
ztype2.add('zt2_P2', value=1.2, vary=False)
ztype2.add('zt2_P3', value=1.3, vary=False)
ztype2.add('zt2_P4', value=1.4, vary=False)

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
    outputl = [0 for i in range(17)]
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
    outarray = odeint(mc.phytomftm_extendedoutput, initcond, timedays_model, args=(all_params, p, z))
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
    elif pfn == 1 and zn == 1:
         print('1P1Z')
         all_params = (standardparams)
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

def test():
    out1P1Z = callmodelrun(1,1)
    out2P2Z = callmodelrun(2,2)
    out3P3Z = callmodelrun(3,3)
    out4P4Z = callmodelrun(4,4)
    out5P5Z = callmodelrun(5,5)

    # plt.subplots_adjust(hspace=0.01)
    f1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 5, sharex='col', sharey='row')

    plotoutput(out1P1Z,1,1,0,'1P1Z')
    plotoutput(out2P1Z,2,2,1,'2P1Z')
    plotoutput(out2P2Z,3,3,2,'2P2Z')
    plotoutput(out2P1Z,4,4,3,'2P1Z')
    plotoutput(out2P2Z,5,5,4,'2P2Z')

    #f1.set_figheight(15)
    #plt.tight_layout()
    plt.show()

    #f1.savefig("foo2.pdf", bbox_inches='tight')

    timedays_model = np.arange(0., 5 * 365., 1.0)

    out1P1Z = callmodelrun(1,1)
    out2P1Z = callmodelrun(2,1)
    out2P3Z = callmodelrun(2,3)
    #out4P4Z = callmodelrun(4,4)

    #out5P5Z = callmodelrun(5,5)
    out2P2Z = callmodelrun(2,2)

    f1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 4, sharex='col', sharey='row')

    plotoutput(out1P1Z,1,1,0, '1P1Z')
    plotoutput(out2P1Z,2,1,1, '2P1Z')
    plotoutput(out2P2Z,2,2,2, '2P2Z')
    plotoutput(out2P3Z,2,3,3, '2P3Z')

timedays_model = np.arange(0., 5 * 365., 1.0)

out4P2Z = callmodelrun(4,2)

f1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 2, sharex='col', sharey='row')

plotoutput(out4P2Z, 4, 2, 1, '4P2Z')