#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
# ExtractEnvFor:
import scipy.interpolate as intrp


# For ParDict:

def sliceparams(pardict, parprefix):
    """Function to extract functioanl type parameters from Parameters object,
    using prefix in key"""
    return {k: v.value for k, v in pardict.items() if k.startswith(parprefix)}


def sliceoffparams(pardict, parprefix):
    """Function to remove e.g. functioanl type parameters from Parameters object,
    using prefix in key"""
    return {k: v.value for k, v in pardict.items() if not k.startswith(parprefix)}


def extractparam(pardict, parsuffix):
    """Function to extract certain single parameter from sliced (!) Parameters object,
    using final characters of key"""
    return next(v for k, v in pardict.items() if k.endswith(parsuffix))


def checkreplaceparam(stdpardict, functypepardict, parsuffix):
    """Function to check that certain single parameter from sliced (!) Parameters object,
    using final characters of key"""
    try:
        ftpara = next(v for k, v in functypepardict.items() if k.endswith(parsuffix))
        return ftpara
    except StopIteration:
        return next(v for k, v in stdpardict.items() if k.endswith(parsuffix))



# For Grazing:

# FUNCTIONS to handle multiple type grazing
# (perhaps implement this within classes, for now outside)

def feedingmatrix(P, Z, pn, zn):
    l = [[phy, zoo] for phy in range(pn) for zoo in range(zn)]
    return l


def grazing(Ptot, Z, pn, zn, zclass):
    l = [[phy, zoo] for phy in range(pn) for zoo in range(zn)]
    all_graze = [[phy, zoo, zclass[zoo].grazing(Ptot, Z[zoo])] for phy, zoo in l]
    return all_graze


def zoogrowth_zoo(GrzingMat, P, pn, zn):
    pairwise = [graz * P[phy] for phy, zoo, graz in GrzingMat]

    sumperzootype = [sum([pairwise[zoo + (phy * zn)] for phy in range(pn)]) for zoo in range(zn)]
    return sumperzootype  # values in list


def zoograzeloss_phy(GrzingMat, pn, zn):
    pairwise = [graz for phy, zoo, graz in GrzingMat]
    sumperphytype = [sum([pairwise[zoo + (phy * zn)] for zoo in range(zn)]) for phy in range(pn)]
    return sumperphytype  # values in list




# FOR FORCING:

# functions adapted from PhytoSFDM Model by Esteban Acevedo-Trejos
def firstderivspl(Forcing, time, k=3, s=None):
    """
    Method to calculate the first derivative of an interpolated spline.

    Parameters
    ----------
    Forcing:
    time: in days
    kind: the type of interpolation either linear, cubic, spline or
           piecewise polynomial
    k: Degree of the smoothing spline
    s: Positive smoothing factor used to choose the number of knots

    Returns
    -------
    The first derivative of the temporally interpolated environmental forcing spline.
    """
    outForcing = Forcing  # spatialave(Forcing)
    tmonth = np.linspace(0., 365., 13) #HERE deprecation warning due to 13. <- float, should be int
    newt = np.mod(time, 365.)
    outintp = intrp.UnivariateSpline(tmonth, outForcing, k=k, s=s)
    return outintp.derivative()(newt)

def dailyinterp(Forcing, time, kind='spline', k=3, s=None):
    """
    Method to interpolate from monthly to daily environmental data.

    Parameters
    -----
    time: in days
    kind: the type of interpolation either linear, cubic, spline or piecewise polynomial
    k: Degree of the smoothing spline
    s: Positive smoothing factor used to choose the number of knots

    Returns
    -------
    The temporally interpolated environmental forcing.
    """
    outForcing = Forcing  # spatialave(Forcing)

    tmonth = np.linspace(0., 12., 13) #HERE again, deprecation warning
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

