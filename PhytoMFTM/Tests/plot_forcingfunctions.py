#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import pandas

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D

from Tests.MODELfunctions_plotParams import all_params, all_params21, all_params22
from PhytoMFTM.ModelClasses import Plankton, Forcing
from lmfit import Parameters, Parameter

from scipy import interpolate




FX = Forcing('constantMLD')


def interpol_deriv():
    Time =  np.linspace(0., 600., 12*50) # np.linspace(0., 1200., 1500)
    inter = FX.MLD.return_interpvalattime(Time)
    deriv = FX.MLD.return_derivattime(Time)
    #deriv2 = FX.MLD.return_derivattime2(Time)
    # Legend
    fig1, (ax1,ax2) = plt.subplots(2)
    ax1.plot(Time, inter, label='inter')
    ax2.plot(Time, deriv, label='deriv')
   # ax3.plot(Time, deriv2, label='deriv')
    #ax4.plot(Time2, full)
    #plt.margins(x=0)
    ax1.set_xlabel('')
    ax1.invert_yaxis()
    ax1.set_ylabel('')
    ax1.set_title('')
    fig1.tight_layout()
    #plt.savefig('PIcurves.png')
    fig1.show()


interpol_deriv()
