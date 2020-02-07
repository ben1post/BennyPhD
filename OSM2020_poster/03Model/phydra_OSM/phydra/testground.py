import pandas
import numpy as np
import scipy.interpolate as intrp
#from phydra.aux import sliceparams, sliceoffparams, checkreplaceparam

from scipy.io import netcdf
import os


def dailyinterp(file, kind='spline', k=3, s=None):
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
    # tmonth = np.linspace(-10.5, 24.473, 12 * 3)
    dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    dpm = dayspermonth * 3
    dpm_cumsum = np.cumsum(dpm) - np.array(dpm) / 2
    if kind == 'spline':
        outintp = intrp.UnivariateSpline(dpm_cumsum, file, k=k, s=s)
        return outintp
    elif kind == 'PWPoly':
        outintp = intrp.PchipInterpolator(dpm_cumsum, file)
        return outintp
    else:
        raise ('Wrong interpolation type passed to dailyinterp function of IndForcing class')


def readCariaco(varname, niskin=True, boxordep='box', forctype='aggTS'): #self
    """read data and return either aggregated regimes or full time series"""
    if niskin:
        niskin_all = pandas.read_csv('Data/NewestData/BoxVSatDepth_02.csv')
        print(list(niskin_all.columns))

    if boxordep == 'box':
        varname_1 = varname+'_Box'
    elif boxordep == 'depth':
        varname_1 = varname+'_AtDepth'

    if forctype == 'aggTS':
        niskin_val = niskin_all[['date', 'month', 'yday', varname_1]]
        niskin_monthly_mean = niskin_val.groupby('month').median()
        # print(niskin_monthly_mean)
        forcing_oneyear = list(niskin_monthly_mean[varname_1])
        forcing_list = forcing_oneyear * 3

    return forcing_list

def return_interpvalattime(y, time):
    """
    Method to return interpolated value of forcing.

    converts time in days to time in months
    """
    newt = np.mod(time, 365.) + 365  # *12./365.
    return y(newt)

def return_derivattime(self, time):
    """
    Method to return derivative (slope) of interpolated value of forcing.

    converts time in days to time in months, and the resulting derivative from per month to per day
    """
    newt = np.mod(time, 365.) + 365  # * 12. / 365.

    if self.forcingtype == "constantMLD" and self.forcvar == "MLD":
        # return 0 as derivative for constant MLD, remove mathematical inaccuracy
        return self.derivative(newt)  # * 0.03
    else:
        return self.derivative(newt)  # * 0.03

x = readCariaco('Temperature')#, boxordep='depth')

print(x)

x_int = dailyinterp(x)

print(x_int)

print(return_interpvalattime(y=x_int, time=100))