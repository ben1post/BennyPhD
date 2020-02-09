#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

from runs.modelrun_aggTS import timedays, outarray, ms

# TODO:
#  - import necessary functions & objects

numcols = 2
f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, numcols, sharex='col')#, sharey='row')


plt.setp((ax1, ax2, ax3, ax4), xticks=[1,60,120,180,240,300,365])
from matplotlib.ticker import MaxNLocator
for axe in (ax1, ax2, ax3, ax4):
    for i in range(numcols):
        axe[i].get_yaxis().set_major_locator(MaxNLocator(nbins=4))
        axe[i].tick_params(top=True, right=True)

# PLOTTING
timedays_ly = timedays[1:366]
outarray_ly = outarray[1460:1825]

# color vectors
#colors = ['#edc951', '#dddddd', '#00a0b0', '#343436', '#cc2a36']
colors = ['#808080','#d55e00', '#cc79a7', '#0072b2', '#009e73', 'grey']
alphas = [1., 0.8, 0.6, 0.4]
lws = [2, 2.5, 4, 5.5]

ax1[0].set_title('model output')

dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
dpm_cumsum = np.cumsum(dayspermonth) - np.array(dayspermonth)/2 #- 15
#print(timedays_ly)


# Figure 1
N = ms.physics.forcing.verif.NO2NO3
#print(N)
# N
N_Max = np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) + np.max(ms.physics.forcing.NOX.return_interpvalattime(timedays)) * 0.1
ax1[0].scatter(N['yday'], N['NO3_NO2_USF_Box'], label='data')
ax1[0].plot(timedays_ly, outarray_ly[:, 0], c=colors[1], lw=lws[0], alpha=alphas[0], label='Model')
ax1[0].set_ylabel('Nutrients \n' '[µM N]', multialignment='center', fontsize=10)
ax1[0].set_ylim(0, N_Max)
ax1[0].legend(fontsize='x-small')


ChlConv = True
# Phyto
CtoChla = 50  # g/g
MolarMassC = 12.0107
CtoNratioPhyto = 6.625
muMolartoChlaconvfactor = CtoChla / MolarMassC / CtoNratioPhyto  # Chla as mg/m-3 to


chla = ms.physics.forcing.verif.FluorChla

chla2 = ms.physics.forcing.verif.HPLC

ax2[0].scatter(chla['yday'], chla['Chlorophyll_Box'] * muMolartoChlaconvfactor, label='Fluor Chla')
#ax2[0].scatter(chla2['yday'], chla2['Tchla'] * muMolartoChlaconvfactor / 100, label='HPLC Chla')


Pall = outarray_ly[:,1]
P_Max = np.max(Pall) + 0.9 * np.max(Pall)

ax2[0].plot(timedays_ly, Pall, c=colors[4], lw=lws[1], label='Model')
ax2[0].legend(fontsize='x-small')
ax2[0].set_ylabel('Phytoplankton \n' '[µM N]', multialignment='center', fontsize=10)
#ax2[0].set_ylim(0, P_Max)

# Z
Zall = outarray_ly[:,2]
Z_Max = np.max(Zall) + 0.1 * np.max(Zall)

ax3[0].plot(timedays_ly, Zall, c=colors[4], lw=lws[1])
ax3[0].set_ylabel('Zooplankton \n' '[µM N]', multialignment='center', fontsize=9)
ax3[0].tick_params('y', labelsize=10)
ax3[0].set_ylim(0, Z_Max)
#ax4[i_plot].set_title('Zooplankton')



# convert PN in µg/L to µM of Detritus!
molarmassNitrogen = 14.0067
mugperlitertomuMolarPN = 1 / molarmassNitrogen  # g/L -> mol/L -> µM

# TODO: add other verification variables! ensure proper converison..
PN = ms.physics.forcing.verif.PN
#ax4[0].scatter(PN['yday'], PN['PN_ug_L_Box'] /14.0067 /100 , label='PN data')

D_Max = np.max(outarray_ly[:, 3]) + 0.2 * np.max(outarray_ly[:, 3])
# D
ax4[0].plot(timedays_ly, outarray_ly[:, 3], c=colors[1], lw=lws[0], alpha=alphas[0])
ax4[0].set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)
ax4[0].set_ylim(0, D_Max)
ax4[0].set_xlabel('Day in year')
# Legend

## PHYSICS ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
muplot = 1
ax1[muplot].set_title('model forcing')

ax4[muplot].set_xlabel('Day in year')

NOX = ms.physics.forcing.NOX.return_interpvalattime(timedays_ly)
NOXdat = ms.physics.forcing.NOX.forcingfile
#print(NOX)
#print(NOXdat)
ax1[muplot].plot(timedays_ly, NOX, c=colors[5], lw=lws[0], alpha=alphas[0])
ax1[muplot].scatter(dpm_cumsum, NOXdat[0:12], c=colors[5])
ax1[muplot].set_ylabel('$N_0$ \n' '[µM]', multialignment='center', fontsize=10)
ax1[muplot].set_ylim(0., N_Max)
#ax1[muplot].invert_yaxis()

MLD = ms.physics.forcing.X258.return_interpvalattime(timedays_ly)
MLDdat = ms.physics.forcing.X258.forcingfile
MLD_max = np.max(MLDdat) + 0.1 * np.max(MLDdat)
ax2[muplot].plot(timedays_ly, MLD, c=colors[5], lw=lws[0], alpha=alphas[0])
ax2[muplot].scatter(dpm_cumsum, MLDdat[0:12], c=colors[5])
ax2[muplot].set_ylabel('$25.8_{iso}$ \n' '[m]', multialignment='center', fontsize=10)
ax2[muplot].set_ylim(0, MLD_max) # 400 for biotrans, 100 for Papa
ax2[muplot].invert_yaxis()

PAR = ms.physics.forcing.PAR.return_interpvalattime(timedays_ly)
PARdat = ms.physics.forcing.PAR.forcingfile
PAR_max = np.max(PAR) + 0.1 * np.max(PAR)
ax3[muplot].plot(timedays_ly, PAR, c=colors[5], lw=lws[0], alpha=alphas[0])
ax3[muplot].scatter(dpm_cumsum, PARdat[0:12], c=colors[5])
ax3[muplot].set_ylabel('PAR \n' '[E $m^{−2}$ $s^{−1}$]', multialignment='center', fontsize=10)
ax3[muplot].set_ylim(0, PAR_max)
# ax1[muplot].invert_yaxis()

Tmld = ms.physics.forcing.SST.return_interpvalattime(timedays_ly)
Tmlddat = ms.physics.forcing.SST.forcingfile
Tmld_max = np.max(Tmld) + 0.1 * np.max(Tmld)
ax4[muplot].plot(timedays_ly, Tmld, c=colors[5], lw=lws[0], alpha=alphas[0])
ax4[muplot].scatter(dpm_cumsum, Tmlddat[0:12], c=colors[5])
ax4[muplot].set_ylabel('$T_{MLD}$ \n' '[°C]', multialignment='center', fontsize=10)
#ax4[muplot].set_ylim(0, Tmld_max)
# ax1[muplot].invert_yaxis()

# Defining custom 'xlim' and 'ylim' values.
xlim = (0, 365)

# Setting the values for all axes.
plt.setp((ax1, ax2, ax3, ax4), xlim=xlim)

f1.align_ylabels()

plt.subplots_adjust(hspace=0.1)

plt.tight_layout()
plt.show()