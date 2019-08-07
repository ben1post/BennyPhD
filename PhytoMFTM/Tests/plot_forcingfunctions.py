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

import scipy.interpolate as intrp


def dailyinterp(outforcing, time, kind='spline', k=5, s=None):
    """
    Method to interpolate from monthly to daily environmental data.

    Parameters
    ----------
    time: in days
    kind: the type of interpolation either linear, cubic, spline or
           piecewise polynomial
    k: Degree of the smoothing spline
    s: Positive smoothing factor used to choose the number of knots

    Returns
    -------
    The temporally interpolated environmental forcing.
    """

    outForcing = np.append(outforcing, outforcing[0])

    tmonth = np.linspace(0., 12., 13.)
    newt = np.mod(time, 365.) * 12. / 365.
    if kind == 'spline':
        outintp = intrp.UnivariateSpline(tmonth, outForcing, k=k, s=s)
        return outintp(newt)
    elif kind == 'PWPoly':
        outintp = intrp.PchipInterpolator(tmonth, outForcing)
        return outintp(newt)
    else:
        outintp = intrp.interp1d(tmonth, outForcing, kind=kind)
        return outintp(newt)


def firstderivspl(outforcing, time, k=5, s=None):
    """
    Method to calculate the first derivative of an interpolated spline.

    Parameters
    ----------
    time: in days
    kind: the type of interpolation either linear, cubic, spline or
           piecewise polynomial
    k: Degree of the smoothing spline
    s: Positive smoothing factor used to choose the number of knots

    Returns
    -------
    The first derivative of the temporally interpolated environmental forcing spline.
    """

    outForcing = np.append(outforcing, outforcing[0])

    tmonth = np.linspace(0., 365., 13.)
    newt = np.mod(time, 365.)
    outintp = intrp.UnivariateSpline(tmonth, outForcing, k=k, s=s)
    return outintp.derivative()(newt)


FX = Forcing('constantMLD')

def interpol_deriv2():
    Time =  np.linspace(0., 600., 12*50) # np.linspace(0., 1200., 1500)
    Time2 =  np.linspace(360., 370., 12*50) # np.linspace(0., 1200., 1500)
    inter = dailyinterp(FX.MLD.forcingfile[0:12], Time)
    inter2 = dailyinterp(FX.MLD.forcingfile[0:12], Time2)
    deriv = firstderivspl(FX.MLD.forcingfile[0:12], Time)
    deriv2 = firstderivspl(FX.MLD.forcingfile[0:12], Time2)
    #deriv2 = FX.MLD.return_derivattime2(Time)
    # Legend
    fig1, (ax1,ax2,ax3,ax4) = plt.subplots(4)
    ax1.plot(Time, inter, label='inter')
    ax1.set_title('Interpolated MLD [m]')
    ax1.set_ylim(0,40)
    ax2.plot(Time2, inter2, label='inter')
    ax2.set_title('Interpolated MLD [m] ZOOM')
    ax2.set_ylim(36, 38)
    ax3.plot(Time, deriv, label='deriv')
    ax3.set_title('Derivative of MLD [m $d^{-1}$]')
    ax3.set_ylim(-0.5, 0.5)
    ax4.plot(Time2, deriv2, label='deriv')
    ax4.set_title('Derivative of MLD [m $d^{-1}$] ZOOM')
    ax4.set_ylim(0., 0.5)
    #ax3.plot(Time, deriv2, label='deriv')
    #ax4.plot(Time2, full)
    #plt.margins(x=0)
    ax1.set_xlabel('')
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    #ax1.set_ylabel('')
    #ax1.set_title('')
    fig1.tight_layout()
    #plt.savefig('PIcurves.png')
    fig1.show()

def interpol_deriv():
    Time =  np.linspace(0., 600., 12*50) # np.linspace(0., 1200., 1500)
    Time2 =  np.linspace(360., 370., 12*50) # np.linspace(0., 1200., 1500)
    inter = FX.MLD.return_interpvalattime(Time)
    inter2 = FX.MLD.return_interpvalattime(Time2)
    deriv = FX.MLD.return_derivattime(Time)
    deriv2 = FX.MLD.return_derivattime(Time2)
    #deriv2 = FX.MLD.return_derivattime2(Time)
    # Legend
    fig1, (ax1,ax2,ax3,ax4) = plt.subplots(4)
    ax1.plot(Time, inter, label='inter')
    ax1.set_title('Interpolated MLD [m]')
    ax1.set_ylim(0,40)
    ax2.plot(Time2, inter2, label='inter')
    ax2.set_title('Interpolated MLD [m] ZOOM')
    ax2.set_ylim(24, 26)
    ax3.plot(Time, deriv, label='deriv')
    ax3.set_title('Derivative of MLD [m $d^{-1}$]')
    ax3.set_ylim(-0.3, 0.4)
    ax4.plot(Time2, deriv2, label='deriv')
    ax4.set_title('Derivative of MLD [m $d^{-1}$] ZOOM')
    ax4.set_ylim(0., 0.5)
    #ax3.plot(Time, deriv2, label='deriv')
    #ax4.plot(Time2, full)
    #plt.margins(x=0)
    ax1.set_xlabel('')
    ax1.invert_yaxis()
    ax2.invert_yaxis()
    #ax1.set_ylabel('')
    #ax1.set_title('')
    fig1.tight_layout()
    #plt.savefig('PIcurves.png')
    fig1.show()


interpol_deriv()
interpol_deriv2()
