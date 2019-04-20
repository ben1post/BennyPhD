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


from Tests.fittingexercise01 import result1,result2,result3,timedays_model

initialcond1 = result1[1]
initialcond2 = result2[1]
initialcond3 = result3[1]

# read yearly data (for comparison to model) from Cariaco
ChlA = pandas.read_csv('../CARIACOModelingProject/Data/DataFiles_Processed/ChlA_bottle_yearly_surface.csv')

NO3NO2 = pandas.read_csv('../CARIACOModelingProject/Data/DataFiles_Processed/NO3NO2_bottle_yearly_surface.csv')

SiOH_UDO = pandas.read_csv('../CARIACOModelingProject/Data/DataFiles_Processed/SiOH_UDO_bottle_yearly_surface.csv')
SiOH_USF = pandas.read_csv('../CARIACOModelingProject/Data/DataFiles_Processed/SiOH_USF_bottle_yearly_surface.csv')

# Zooplankton:
ZooBM = pandas.read_csv('../CARIACOModelingProject/Data/DataFiles_Processed/ZooBM_All.csv')
# zooplankton biomass is in [mg/m^3 dry weight]
# for now, a 10% N content of the dry weight of all Zooplankton is assumed,
# [Gorsky et al. 1988] "C and N composition of some northwestern Mediterranean zooplankton and micronekton species"
# to convert the values to µM  N  ::  1 mg N/m^3 = 0.071394 μM N

###########--monthly medians---######
ChlA = ChlA.assign(month = pandas.to_datetime(ChlA['yday'], format='%j').dt.month)
NO3NO2 = NO3NO2.assign(month = pandas.to_datetime(NO3NO2['yday'], format='%j').dt.month)
SiOH_USF = SiOH_USF.assign(month = pandas.to_datetime(SiOH_USF['yday'], format='%j').dt.month)
ZooBM = ZooBM.assign(month = pandas.to_datetime(ZooBM['yday'], format='%j').dt.month)

ChlA_monthly_median = ChlA.groupby('month').median()
NO3NO2_monthly_median = NO3NO2.groupby('month').median()
SiOH_USF_monthly_median = SiOH_USF.groupby('month').median()
ZooBM_monthly_median = ZooBM.groupby('month').median()



# color vectors
colors = ['#00C90D', '#01939A', '#d16f00', '#d13700', '#d10000']
alphas = [1., 0.8, 0.6, 0.4]
lws = [1, 2.5, 4, 5.5]

# artist for legends
FullArtist = plt.Line2D((0, 1), (0, 0), c=colors[4], alpha=alphas[1], lw=lws[0])

def plot(pfn,zn,result):
    # Figure 1

    outarray = g(initialcond1, timedays_model, result.params)

    timedays = timedays_model[1:366]
    # truncate outarraySiNO to last year of 5:
    outarray_ly = outarray[1460:1825]


    f1, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(5, 1, sharex='col', sharey='row')
    # N
    ax1.plot(timedays, outarray_ly[:, 0], c=colors[4], lw=lws[0], alpha=alphas[0], label='Model')
    ax1.set_ylabel('Nitrate \n' '[µM]', multialignment='center', fontsize=10)
    ax1.set_ylim(-0.1, 5)

    # N Data
    ax1.scatter(NO3NO2['yday'].values, NO3NO2['NO3NO2'].values, c=colors[1],s=4.3, label='Data')
    ax1.set_title('Nitrate & NO3 + NO2')
    ax1.legend(loc=1)

    # Si
    ax2.plot(timedays, outarray_ly[:, 1], c=colors[4], lw=lws[0], alpha=alphas[0])
    ax2.set_ylabel('Silicate \n' '[µM]', multialignment='center', fontsize=10)
    ax2.set_ylim(-0.1, 12)
    # Si Data
    ax2.scatter(SiOH_USF['yday'].values, SiOH_USF['SiOH'].values, c=colors[1],s=4.3)
    ax2.set_title('Silicate & SiOH')

    #Phyto
    [ax3.plot(timedays, outarray_ly[:, 3+zn+i], c=colors[i+2]) for i in range(pfn)]
    ax3.set_ylabel('Phyto \n' '[µM N]', multialignment='center', fontsize=10)
    ax3.set_ylim(-0.1, 0.8)
    # Phyto Data
    ax3_tx = ax3.twinx()
    ax3_tx.scatter(ChlA['yday'].values, ChlA['ChlA'].values, c=colors[1],s=4.3)
    ax3_tx.set_ylabel('ChlA \n [mg/m3]')

    ax3.set_title('Phy Biomass & ChlA Data')

    # Z
    [ax4.plot(timedays, outarray_ly[:, 3+i], c=colors[i+3], lw=lws[0], alpha=alphas[0]) for i in range(pfn)]
    ax4.set_ylabel('Zooplankton \n' '[µM N]', multialignment='center', fontsize=9)
    ax4.tick_params('y', labelsize=10)

    ax4.scatter(ZooBM['yday'].values, ZooBM['ZooBM'].values*0.1*0.071394, c=colors[1],s=4.3)
    # 10% N of dry weight assumed, converted to µM

    ax4.set_title('Zooplankton')
    ax4.set_ylim(0, 0.62)

    # D
    ax5.plot(timedays, outarray_ly[:, 3], c=colors[1], lw=lws[0], alpha=alphas[0])
    ax5.set_ylabel('Detritus \n' '[µM N]', multialignment='center', fontsize=9)

    ax5.set_title('Detritus')
    ax5.set_ylim(bottom=0)

    ax5.set_xlabel('Day in year', fontsize=14)
    # Legend


    #plt.subplots_adjust(hspace=0.01)
    f1.set_figheight(15)
    plt.tight_layout()
    plt.show()


plot(1,1,result1)
plt.show()