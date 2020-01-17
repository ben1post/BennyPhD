#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import pandas
import matplotlib.pyplot as plt
from PhytoMFTM.ModelClasses import Forcing


def returnintverifdat(DF):
    DF = DF.assign(month=pandas.to_datetime(DF['yday'], format='%j').dt.month)
    DF_monthly_median = DF.groupby('month').median()

    DF_m1 = DF_monthly_median.copy()
    DF_m2 = DF_monthly_median.copy()
    DF_m3 = DF_monthly_median.copy()

    DF_m2['yday'] = DF_m2['yday'] + 365
    DF_m3['yday'] = DF_m3['yday'] + 2 * 365

    DF_x = pandas.concat([DF_m1, DF_m2, DF_m3])

    tm_dat_conc = np.arange(0., 3 * 365., 1.0)

    DF_pad = DF_monthly_median.set_index('yday').reindex(tm_dat_conc).reset_index()

    DF_int = DF_pad.Value.interpolate().values  # .plot() # .ffill().bfill()
    DF_intslice = DF_int[365:365 + 365]

    tm_year = np.arange(0., 365., 1.0)

    DF_int2 = pandas.DataFrame()
    DF_int2['yday'] = tm_year
    DF_int2['Value'] = DF_intslice
    return DF_int2


fx = Forcing('constantMLD')

print(fx)

modeldepth = 100
CtoNratioPhyto = 6.625
CtoChla = 50
# model output is: µM N - convert to ng m^⁻2 / µM is µmol/l, *1000 = mol/l
muMolartoChlaconvfactor = 1000 * 1000 * modeldepth * CtoNratioPhyto / CtoChla * 891 / 1000000000  # * 90 # µM N m^-2 * C N^-1 / C chla^-1 = µM chla m^-2 to ng m^-2

Tchla = pandas.read_csv('ValidationData/Tchla_r1.csv')
Tchla2 = pandas.read_csv('ValidationData/Tchla_r2.csv')

Tchla['Date'] = pandas.to_datetime(Tchla.Date, format='%Y-%m-%d')
Tchla2['Date'] = pandas.to_datetime(Tchla2.Date, format='%Y-%m-%d')

Tchla = Tchla[Tchla['spec']=='Tchla']
Tchla2 = Tchla2[Tchla2['spec']=='Tchla']

Tchla_monthly_median = Tchla.groupby('month').median()

X21_r1 = pandas.read_csv('Forcing/X21Iso/X21Iso_r1.csv')

#print(Tchla)
#Tchla_int = returnintverifdat(Tchla)
#Tchla2_int = returnintverifdat(Tchla2)


timedays_model = np.arange(0., 365., 1.0)

numplots = 2
f1, (ax1, ax2, ax3, ax4) = plt.subplots(4, numplots, sharex='col')

plt.setp((ax1, ax2), xticks=[1, 60, 120, 180, 240, 300, 360])
from matplotlib.ticker import MaxNLocator
for axe in (ax1, ax2):
    for i in range(numplots):
        axe[i].get_yaxis().set_major_locator(MaxNLocator(nbins=4))

int_MLD = np.array(fx.MLD.return_interpvalattime(timedays_model))

#### MLD ####
ax1[0].plot(timedays_model,int_MLD)
ax1[0].set_ylim(bottom=0)
ax1[0].set_title('MLD')

ax1[1].plot(timedays_model,int_MLD)
ax1[1].invert_yaxis()
ax1[1].set_ylim(top=0)
ax1[1].set_title('MLD - inverse')

deriv_MLD = np.array(fx.MLD.return_derivattime(timedays_model))
ax2[0].plot(timedays_model,deriv_MLD)
ax2[0].set_title('MLD deriv')

K = 0.1 + np.max(deriv_MLD, 0) / np.array(int_MLD)  # i.e. there is constant mixing & increased loss with MLD shallowing

K_Z = deriv_MLD / int_MLD  # i.e. concentration varies with the MLD depth

ax2[1].plot(timedays_model,K)
ax2[1].set_ylim(bottom=0)
#ax2[1].set_ylim(0,1)
ax2[1].set_title('K')

##### NUTRIENTS! ######

#NOx = np.array(fx.NOX.return_interpvalattime(timedays_model))
#ax3[0].plot(timedays_model, NOx)
#ax3[0].set_ylim(bottom=0)
#ax3[0].set_title('N')


#ax3[1].plot(Tchla_int['yday'].values, Tchla_int['Value'].values, color='black')
ax3[1].scatter(Tchla['yday'].values, Tchla['Value'].values, color='green',alpha=0.5)
ax3[1].scatter(Tchla_monthly_median['yday'].values,Tchla_monthly_median['Value'].values)
ax3[1].set_ylim(0,180)
ax3[1].set_title('Total Chl a')


ax3[0].scatter(X21_r1['yday'].values, X21_r1['depth'].values)
ax3[0].set_ylim(0,150)
ax3[0].invert_yaxis()
ax3[0].set_title('21°C - Isopleth (RAW)')

X21 = np.array(fx.X21.return_interpvalattime(timedays_model))
ax4[0].plot(timedays_model, X21)
ax4[0].set_ylim(0,150)
ax4[0].invert_yaxis()
ax4[0].set_title('21°C - Isopleth (interpolated)')

def mixing(x21):
    return 0.1/ x21 * 100

K = (0.1 + np.max(deriv_MLD, 0)) / np.array(int_MLD)  # i.e. there is constant mixing & increased loss with MLD shallowing

ax4[1].plot(timedays_model, mixing(X21))
#ax4[1].set_ylim(0,1)
ax4[1].set_ylim(bottom=0)
#ax4[1].invert_yaxis()
#ax4[1].set_title('21°C - Isopleth (interpolated)')

plt.show()



"""
print(np.mean(Tchla['val'].values))
print(np.mean(Tchla2['val'].values))

plt.scatter(Tchla['Date'].values, Tchla['val'].values)
plt.scatter(Tchla2['Date'].values, Tchla2['val'].values)
plt.show()
"""