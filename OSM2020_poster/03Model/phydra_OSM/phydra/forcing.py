#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas
import numpy as np
import scipy.interpolate as intrp
#from phydra.aux import sliceparams, sliceoffparams, checkreplaceparam

from scipy.io import netcdf
import os

# TODO:
#  - read new TS forcing
#  - add (old) aggTS Regimes to run model in both regimes at CARIACO
#  - read new full TS forcing

class VerifData:
    """
    initializes and reads verification data, contained in ncdf files
    """
    def __init__(self, Lat, Lon, RBB):
        self.Lat = Lat
        self.Lon = Lon
        self.RangeBB = RBB
        self.fordir = os.path.split(os.path.realpath(__file__))[0]
        self.chla = self.readchla()
        self.chlaint = self.interplt(self.chla)
        self.N = self.readnutrientaboveMLD('n0')
        self.Nint = self.interplt(self.N)

        print('VerifData forcing created')

    def readchla(self):
        ncfile = netcdf.netcdf_file(self.fordir + '/ChlAclimatology_MODISaqua_L3_nc3.nc', 'r')
        nclat = ncfile.variables['lat'].data.copy()
        nclon = ncfile.variables['lon'].data.copy()
        ncdat = ncfile.variables['chlor_a'].data.copy()
        ncfile.close()
        longrid, latgrid = np.meshgrid(nclon, nclat)
        selectarea = np.logical_and(longrid <= self.Lon + self.RangeBB, longrid >= self.Lon - self.RangeBB) * \
                     np.logical_and(latgrid <= self.Lat + self.RangeBB, latgrid >= self.Lat - self.RangeBB)
        outforcing = list(np.nanmean(ncdat[:, selectarea], axis=1))
        return outforcing

    def readnutrientaboveMLD(self,var):
        if var == 'n0':
            ncfile = netcdf.netcdf_file(self.fordir + '/NitrateAboveMLD_WOA_tmld_test01_nc3.nc', 'r')
        elif var == 'p0':
            ncfile = netcdf.netcdf_file(self.fordir + '/Phosphate_WOA_tmld_test03_nc3.nc', 'r')
        elif var == 'si0':
            ncfile = netcdf.netcdf_file(self.fordir + '/SilicateAboveMLD_WOA_tmld_test02_nc3.nc', 'r')

        nclat = ncfile.variables['lat'].data.copy()
        nclon = ncfile.variables['lon'].data.copy()
        ncdat = ncfile.variables[var].data.copy()
        ncfile.close()
        longrid, latgrid = np.meshgrid(nclon, nclat)
        selectarea = np.logical_and(latgrid <= self.Lat + self.RangeBB, latgrid >= self.Lat - self.RangeBB) * \
                     np.logical_and(longrid <= self.Lon + self.RangeBB, longrid >= self.Lon - self.RangeBB)
        outforcing = list(np.nanmean(ncdat[selectarea, :], axis=0))
        return outforcing

    def interplt(self,dat):
        dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        dpm_3 = dayspermonth * 3
        dpm_cumsum = np.cumsum(dpm_3) - np.array(dpm_3)/2
        dat_dpm = pandas.DataFrame(np.column_stack([dat*3, dpm_cumsum]), columns=['Value', 'yday'])
        tm_dat_conc = np.arange(0., 3 * 365., 1.0)

        dat_pad = dat_dpm.set_index('yday').reindex(tm_dat_conc).reset_index()
        dat_int = dat_pad.Value.interpolate().values[365:365+365]
        return dat_int


class CARIACOdata:
    """
    initializes and reads new forcing & verification data of full time series & supplies this to model
    either aggregated per regime or as full time series
    """
    def __init__(self, forcvar, data='niskin', time='regime1', k=3, s=None, kind="spline", boxordep='box', forctype='aggTS'):
        """# parameters for interpolation"""
        self.k = k
        self.s = s
        self.kind = kind
        self.forcingtype = forctype
        self.forcvar = forcvar

        self.forcingfile = self.readCariaco(forcvar, data=data, time=time, boxordep=boxordep, forctype=forctype)

        self.interpolated = self.dailyinterp(self.forcingfile, self.kind, self.k, self.s)

        if kind == "spline":
            self.derivative = self.interpolated.derivative()
            self.derivative2 = self.interpolated.derivative(n=2)
        print(forcvar + ' forcing created')

    def readCariaco(self, varname, data, time, boxordep, forctype):  # self
        """read data and return either aggregated regimes or full time series"""
        if data == 'niskin':
            df_all = pandas.read_csv('Data/NewestData/BoxVSatDepth_02.csv')
            if boxordep == 'box': varname = varname + '_Box'
            elif boxordep == 'depth': varname = varname + '_AtDepth'

        elif data == 'SeaWiFS': df_all = pandas.read_csv('Data/NewestData/PARXSeaWiFS_02.csv')

        elif data == 'x25.8': df_all = pandas.read_csv('Data/NewestData/x258_02.csv')

        if forctype == 'aggTS':
            df_val = df_all[['date', 'month', 'yday', varname]]
            if time == 'regime1':
                df_series = df_val.set_index(df_val['date'])
                df_series_cut = df_series.loc['1996-01-01':'2000-10-30']
                df_val = df_series_cut
                print(df_val)
            elif time == 'regime2':
                df_series = df_val.set_index(df_val['date'])
                df_series_cut = df_series.loc['2006-06-30':'2010-12-31']
                df_val = df_series_cut
                print(df_val)
            df_monthly_mean = df_val.groupby('month').mean()
            # print(niskin_monthly_mean)
            forcing_oneyear = list(df_monthly_mean[varname])
            forcing_list = forcing_oneyear * 3

        return forcing_list

    def dailyinterp(self, file, kind, k, s):
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

    def return_interpvalattime(self, time):
        """
        Method to return interpolated value of forcing.

        converts time in days to time in months
        """
        newt = np.mod(time, 365.) + 365  # *12./365.
        return self.interpolated(newt)

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


class Forcing:
    """
    Class to initialze all other forcings, and read files,
    call interpolation on subclasses
    """
    def __init__(self, forcingtype, time):
        if forcingtype == 'aggTS':
            self.X258 = CARIACOdata('x258depth', data='x25.8', boxordep='depth', time=time, k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = CARIACOdata('NO3_NO2_USF', data='niskin', time=time, k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = CARIACOdata('Temperature', data='niskin', time=time, k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = CARIACOdata('value', data='SeaWiFS', time=time, k=5, s=None, kind="spline", forctype=forcingtype)
            self.type = 'Box'

        else:
            raise('wrong forcingtype passed to Forcing class')