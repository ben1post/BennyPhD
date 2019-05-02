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
standardparams.add('pfun_num', value=1, vary=False) #, min=1,max=2,brute_step=1)
# number of zooplankton groups
standardparams.add('zoo_num', value=1, vary=False) #, min=1,max=2,brute_step=1)
# mld - related
standardparams.add('kappa', value=0.1, vary=False) #min=0.09, max=0.11)      # Diffusive mixing across thermocline (m*d^-1)
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
"""
add feeding pref params!
Z1P1
Z1P2
Z2P1
Z2P2
"""

# set up phytoplankton type 1 (e.g. diatoms)
ptype1 = Parameters()
ptype1.add('pt1_U_N', value=5.0, vary=False)    # , min=1, max=6)#, vary=False)   # Silicate Half Saturation Constant
ptype1.add('pt1_muP', value=1.4, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype1.add('pt1_ratioSi', value=0., vary=False)  # Silicate ratio


# set up phytoplankton type 2 (e.g. other)
ptype2 = Parameters()
ptype2.add('pt2_U_N', value=4.3, vary=False)     # Silicate Half Saturation Constant
ptype2.add('pt2_muP', value=1.1, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype2.add('pt2_ratioSi', value=0., vary=False)  # Silicate ratio

# set up phytoplankton type 3 (e.g. other)
ptype3 = Parameters()
ptype3.add('pt3_U_Si', value=0., vary=False)     # Silicate Half Saturation Constant
ptype3.add('pt3_U_N', value=1.94, vary=False)    # Nitrate Half Saturation Constant
ptype3.add('pt3_muP', value=0.5, vary=False)    # Phytoplankton maximum growth rate (d^-1)
ptype3.add('pt3_ratioSi', value=0., vary=False)  # Silicate ratio



# set up zooplankton type 1 (e.g. small zooplankton)
ztype1 = Parameters()
ztype1.add('zt1_gr_p', value=0.6, vary=False)    # , min=0.4, max=1.)#False)   # Portion of Phytoplankton being grazed by Zooplankton
ztype1.add('zt1_muZ', value=0.1, vary=False)    # Zooplankton maximum grazing rate (d^-1)

#ztype1.add('zt1_moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
#ztype1.add('zt1_deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)


# set up zooplankton type 2 (e.g. larger zooplankton)
ztype2 = Parameters()
ztype2.add('zt2_gr_p', value=0.6, vary=False)   # Portion of Phytoplankton being grazed by Zooplankton
ztype2.add('zt2_muZ', value=0.1, vary=False)    # Zooplankton maximum grazing rate (d^-1)

#ztype2.add('zt2_moZ', value=0.1, vary=False)        # Zooplankton mortality (d^-1)
#ztype2.add('zt2_deltaZ', value=0.31, vary=False)    # Zooplankton Grazing assimilation coefficient (-)


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

    # number of phytoplankton func types
    pfn = params['pfun_num'].value
    # number of zooplankton groups
    zn = params['zoo_num'].value

    if pfn == 2 and zn == 2:
        print('2P2Z')
        all_params = (standardparams + ptype1 + ptype2 + ztype1 + ztype2)
    elif pfn == 2 and zn == 1:
        print('2P1Z')
        all_params = (standardparams + ptype1 + ptype2)
    elif pfn == 1 and zn == 1:
        print('1P1Z')
        all_params = (standardparams)
    else:
        print('just standard params')
        all_params = (standardparams)

    all_params = (standardparams)
    parameters = all_params
    initialcond = setupinitcond(pfn, zn)
    print(initialcond)

    z = Plankton(all_params, 'Zooplankton').init()
    p = Plankton(all_params, 'Phytoplankton').init()

    # INTEGRATE:
    tos = time.time()
    print('starting integration')
    outarray = odeint(mc.phytomftm_extendedoutput, initcond, timedays_model, args=(all_params, p, z))
    tos1 = time.time()
    print('finished after %4.3f sec' % (tos1 - tos))

    return outarray


def callmodelrun(params):
    # number of phytoplankton func types
    pfn = params['pfun_num'].value
    # number of zooplankton groups
    zn = params['zoo_num'].value

    if pfn == 2 and zn == 2:
         print('2P2Z')
         all_params = (standardparams + ptype1 + ptype2 + ztype1 + ztype2)
    elif pfn == 2 and zn == 1:
         print('2P1Z')
         all_params = (standardparams + ptype1 + ptype2)
    elif pfn == 1 and zn == 1:
         print('1P1Z')
         all_params = (standardparams)
    else:
        print('just standard params')
        all_params = (standardparams)

    all_params = (standardparams)
    parameters = all_params
    initialcond = setupinitcond(pfn,zn)
    print(initialcond)
    out = runmodel(parameters,initialcond)

    return out



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




def residual(paras, *initialcond):
    """
    compute the residual between actual data and fitted data
    """
    pfn = paras['pfun_num'].value
    zn = paras['zoo_num'].value

    model = g(initialcond, timedays_model, paras)


    # to implement fitting algorithm make sure to calculate residual only for the last year!

    # will have to 1. : simplify the data (i.e. median per month)
    # will have to put data into structure to calculate efficiently (i.e. pandas dataframe like df[1] = N, df[2] = Si, etc.)
    model_ly = model[1460:1825]

    # aggregate model output in the same way as validation data (monthly mean)
    # create month vector to add to model output dataframe for analysis
    oneyearmodel = pandas.DataFrame()
    oneyearmodel = oneyearmodel.assign(day=pandas.Series(np.linspace(1, 365, 365)))

    # combine two columns
    phyto_model = pandas.DataFrame(
        {'data':sum([model_ly[:,3+zn+i] for i in range(pfn)]), 'month': pandas.to_datetime(oneyearmodel['day'], format='%j').dt.month})
    phyto_monthly_median = phyto_model.groupby('month').median()
    phyto_resid = (phyto_monthly_median['data'].values - ChlA_monthly_median['ChlA'].values * 0.1)

    nitrate_model = pandas.DataFrame(
        {'data': model_ly[:, 0], 'month': pandas.to_datetime(oneyearmodel['day'], format='%j').dt.month})
    nitrate_monthly_median = nitrate_model.groupby('month').median()
    nitrate_resid = (nitrate_monthly_median['data'].values - NO3NO2_monthly_median['NO3NO2'].values * 0.1)

    silicate_model = pandas.DataFrame(
        {'data': model_ly[:, 1], 'month': pandas.to_datetime(oneyearmodel['day'], format='%j').dt.month})
    silicate_monthly_median = silicate_model.groupby('month').median()
    silicate_resid = (silicate_monthly_median['data'].values - SiOH_USF_monthly_median['SiOH'].values * 0.1)

    zoo_model = pandas.DataFrame(
        {'data': sum([model_ly[:, 3 + i] for i in range(zn)]), 'month': pandas.to_datetime(oneyearmodel['day'], format='%j').dt.month})
    zoo_monthly_median = zoo_model.groupby('month').median()
    zoo_resid = (zoo_monthly_median['data'].values - ZooBM_monthly_median['ZooBM'].values * 0.1)

    ss = np.concatenate((phyto_resid, nitrate_resid, silicate_resid, zoo_resid))
    return ss

def g(x0, t, paras):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    z = Plankton(paras, 'Zooplankton').init()
    p = Plankton(paras, 'Phytoplankton').init()

    x = odeint(mc.phytomftm_extendedoutput, x0, t, args=(paras,p,z))
    return x


def do_fit(pfn,zn):
    standardparams.add('pfun_num', value=pfn, vary=False)  # , min=1,max=2,brute_step=1)
    # number of zooplankton groups
    standardparams.add('zoo_num', value=zn, vary=False)  # , min=1,max=2,brute_step=1)
    if pfn == 2 and zn == 2:
        print('2P2Z')
        all_params = (standardparams + ptype1 + ptype2 + ztype1 + ztype2)
    elif pfn == 2 and zn == 1:
        print('2P1Z')
        all_params = (standardparams + ptype1 + ptype2)
    elif pfn == 1 and zn == 1:
        print('1P1Z')
        all_params = (standardparams)
    else:
        print('just standard params')
        all_params = (standardparams)

    #all_params = (standardparams)
    print(pfn,zn)
    initialcond = setupinitcond(pfn, zn)

    print(initialcond)

    result = minimize(residual, all_params, args=initialcond, method='least_sq')  # leastsq nelder
    # check results of the fit
    #outarray = g(initcond, timedays_model, result.params)
    print(result.aic)

    report_fit(result)
    return result, initialcond


timedays_model = np.arange(0., 5 * 365., 1.0)


result1 = do_fit(1,1)
# result2 = do_fit(2,1)
# result3 = do_fit(2,2)