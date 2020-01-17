#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas
import numpy as np
import scipy.interpolate as intrp
from PhytoMFTM.AuxFuncs import sliceparams, sliceoffparams, checkreplaceparam

from scipy.io import netcdf
import os

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


class WOAForcing:
    """
    initializes and reads forcing from a certain location in the WOA 2009 data, contained in ncdf files

    """
    def __init__(self, lat, lon, rangebb, varname):
        self.Lat = lat
        self.Lon = lon
        self.RangeBB = rangebb
        self.varname = varname
        self.fordir = os.path.split(os.path.realpath(__file__))[0]
        self.outForcing = self.spatialave()

    def spatialave(self):
        """
        Method to extract spatially averaged environmental forcing.

        Returns
        -------
        The spatial average of the respective environmental forcing per month.
        """

        if self.varname == 'mld':
            ncfile = netcdf.netcdf_file(self.fordir + '/mld_mindtr02_l3_nc3.nc', 'r')
            # print(ncfile.dimensions)
            nclat = ncfile.variables['lat'].data.copy()
            nclon = ncfile.variables['lon'].data.copy()
            ncdat = ncfile.variables['mld_mindtr02_rmoutliers_smth_okrg'].data.copy()
            ncfile.close()
            longrid, latgrid = np.meshgrid(nclon, nclat)
            selectarea = np.logical_and(longrid <= self.Lon + self.RangeBB, longrid >= self.Lon - self.RangeBB) * \
                         np.logical_and(latgrid <= self.Lat + self.RangeBB, latgrid >= self.Lat - self.RangeBB)
            outforcing = list(np.nanmean(ncdat[:, selectarea], axis=1))
            return outforcing * 3

        elif self.varname == 'par':
            ncfile = netcdf.netcdf_file(self.fordir + '/PARclimatology_MODISaqua_L3_nc3.nc', 'r')
            nclat = ncfile.variables['lat'].data.copy()
            nclon = ncfile.variables['lon'].data.copy()
            ncdat = ncfile.variables['par'].data.copy()
            ncfile.close()
            longrid, latgrid = np.meshgrid(nclon, nclat)
            selectarea = np.logical_and(longrid <= self.Lon + self.RangeBB, longrid >= self.Lon - self.RangeBB) * \
                         np.logical_and(latgrid <= self.Lat + self.RangeBB, latgrid >= self.Lat - self.RangeBB)
            outforcing = list(np.nanmean(ncdat[:, selectarea], axis=1))
            return outforcing * 3

        elif self.varname == 'n0x':
            ncfile = netcdf.netcdf_file(self.fordir + '/Nitrate_WOA_tmld_test03_nc3.nc', 'r')
            nclat = ncfile.variables['lat'].data.copy()
            nclon = ncfile.variables['lon'].data.copy()
            ncdat = ncfile.variables['n0'].data.copy()
            ncfile.close()
            longrid, latgrid = np.meshgrid(nclon, nclat)
            selectarea = np.logical_and(latgrid <= self.Lat + self.RangeBB, latgrid >= self.Lat - self.RangeBB) * \
                         np.logical_and(longrid <= self.Lon + self.RangeBB, longrid >= self.Lon - self.RangeBB)
            outforcing = list(np.nanmean(ncdat[selectarea, :], axis=0))
            return outforcing * 3

        elif self.varname == 'p0x':
            ncfile = netcdf.netcdf_file(self.fordir + '/Phosphate_WOA_tmld_test03_nc3.nc', 'r')
            nclat = ncfile.variables['lat'].data.copy()
            nclon = ncfile.variables['lon'].data.copy()
            ncdat = ncfile.variables['p0'].data.copy()
            ncfile.close()
            longrid, latgrid = np.meshgrid(nclon, nclat)
            selectarea = np.logical_and(latgrid <= self.Lat + self.RangeBB, latgrid >= self.Lat - self.RangeBB) * \
                         np.logical_and(longrid <= self.Lon + self.RangeBB, longrid >= self.Lon - self.RangeBB)
            outforcing = list(np.nanmean(ncdat[selectarea, :], axis=0))
            return outforcing * 3

        elif self.varname == 'si0x':
            ncfile = netcdf.netcdf_file(self.fordir + '/Silicate_WOA_tmld_test02_nc3.nc', 'r')
            nclat = ncfile.variables['lat'].data.copy()
            nclon = ncfile.variables['lon'].data.copy()
            ncdat = ncfile.variables['si0'].data.copy()
            ncfile.close()
            longrid, latgrid = np.meshgrid(nclon, nclat)
            selectarea = np.logical_and(latgrid <= self.Lat + self.RangeBB, latgrid >= self.Lat - self.RangeBB) * \
                         np.logical_and(longrid <= self.Lon + self.RangeBB, longrid >= self.Lon - self.RangeBB)
            outforcing = list(np.nanmean(ncdat[selectarea, :], axis=0))
            return outforcing * 3

        elif self.varname == 'sst':
            ncfile = netcdf.netcdf_file(self.fordir + '/Temp_WOA_tmld_test02_nc3.nc', 'r')
            nclat = ncfile.variables['lat'].data.copy()
            nclon = ncfile.variables['lon'].data.copy()
            ncdat = ncfile.variables['t_mld'].data.copy()
            ncfile.close()
            longrid, latgrid = np.meshgrid(nclon, nclat)
            selectarea = np.logical_and(latgrid <= self.Lat + self.RangeBB, latgrid >= self.Lat - self.RangeBB) * \
                         np.logical_and(longrid <= self.Lon + self.RangeBB, longrid >= self.Lon - self.RangeBB)
            outforcing = list(np.nanmean(ncdat[selectarea, :], axis=0))
            return outforcing * 3
        else:
            return 'Please specify either mld, par, n0x or sst'



class IndForcing:
    """
    initializes and reads individual model forcings, and contains methods for daily interpolation and derivation

    """
    def __init__(self, forcvar, filepath, k=3, s=None, kind="spline", forctype=None, WOA=False, Lat=47.5,Lon=-15.5,RBB=2.5):
        # parameters for interpolation
        self.k = k
        self.s = s
        self.kind = kind
        self.forcingtype = forctype
        self.forcvar = forcvar
        if (WOA == False) and (forctype != 'EMPOWER'):
            self.forcingfile = self.readconcforc(forcvar, filepath)
        elif forctype == 'EMPOWER':
            self.forcingfile = self.readEMPOWER(forcvar, filepath)
        else:
            self.forcingfile = WOAForcing(Lat, Lon, RBB, forcvar).outForcing
        self.interpolated = self.dailyinterp(self.forcingfile, self.kind, self.k, self.s)
        if kind == "spline":
            self.derivative = self.interpolated.derivative()
            self.derivative2 = self.interpolated.derivative(n=2)


        print(forcvar + ' forcing created')
    def readEMPOWER(self, varname, filepath):
        forc_all = pandas.read_csv(filepath)
        forcing_monthly_median = forc_all.groupby('month').mean()
        forcing_oneyear = list(forcing_monthly_median[varname])
        forcing_list = forcing_oneyear * 3
        return forcing_list


    def readconcforc(self, varname, filepath):
        """ read forcing from csv file and calculate monthly means """
        forc = pandas.read_csv(filepath)
        forcing_monthly_median = forc.groupby('month').mean()
        forcing_oneyear = list(forcing_monthly_median[varname])
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
        dpm_cumsum = np.cumsum(dpm) - np.array(dpm)/2
        if kind == 'spline':
            outintp = intrp.UnivariateSpline(dpm_cumsum, file, k=k, s=s)
            return outintp
        elif kind == 'PWPoly':
            outintp = intrp.PchipInterpolator(dpm_cumsum, file)
            return outintp
        else:
            raise('Wrong interpolation type passed to dailyinterp function of IndForcing class')

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
    def __init__(self, forcingtype, lat=0, lon=0, rbb=1):
        if forcingtype == 'variableMLD':
            self.MLD = IndForcing('MLD', 'Forcing/MLDdriven/MLD2015_R1.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = IndForcing('NO2NO3', 'Forcing/MLDdriven/NO2NO3_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SiOX = IndForcing('SiOH', 'Forcing/MLDdriven/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/MLDdriven/SST_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/MLDdriven/PAR_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.type = 'MLD'

        elif forcingtype == 'varMLDconstNuts':
            self.MLD = IndForcing('MLD', 'Forcing/MLDdriven/MLD2015_R1.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = IndForcing('NO2NO3', 'Forcing/constantMLD/NO2NO3_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SiOX = IndForcing('SiOH', 'Forcing/constantMLD/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/MLDdriven/SST_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/MLDdriven/PAR_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.type = 'MLD'

        elif forcingtype == 'constantKappa':
            self.MLD = IndForcing('MLD', 'Forcing/MLDdriven/MLD2015_R1.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = IndForcing('NO2NO3', 'Forcing/constantMLD/NO2NO3_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SiOX = IndForcing('SiOH', 'Forcing/constantMLD/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/constantMLD/SST_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/constantMLD/PAR_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.type = 'box_constantKappa'

        elif forcingtype == 'stochasticKappa':
            self.MLD = IndForcing('MLD', 'Forcing/MLDdriven/MLD2015_R1.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = IndForcing('NO2NO3', 'Forcing/constantMLD/NO2NO3_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SiOX = IndForcing('SiOH', 'Forcing/constantMLD/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/constantMLD/SST_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/constantMLD/PAR_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.type = 'box_stochasticKappa'

        elif forcingtype == 'MLD_and_stochastic':
            self.MLD = IndForcing('MLD', 'Forcing/MLDdriven/MLD2015_R1.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = IndForcing('NO2NO3', 'Forcing/constantMLD/NO2NO3_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SiOX = IndForcing('SiOH', 'Forcing/constantMLD/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/constantMLD/SST_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/constantMLD/PAR_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.type = 'box_MLD_stochastic'

        elif forcingtype == 'constantMLD':
            self.MLD = IndForcing('MLD', 'Forcing/MLDdriven/MLD2015_R1.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = IndForcing('NO2NO3', 'Forcing/constantMLD/NO2NO3_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SiOX = IndForcing('SiOH', 'Forcing/constantMLD/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/constantMLD/SST_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/constantMLD/PAR_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.X21 = IndForcing('depth', 'Forcing/X21Iso/X21Iso_r1.csv', k=3, s=1000, kind="spline",
                                  forctype=forcingtype)
            self.type = 'box'

        elif forcingtype == 'flowthrough':
            self.MLD = IndForcing('MLD', 'Forcing/MLDdriven/MLD2015_R1.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = IndForcing('NO2NO3', 'Forcing/constantMLD/NO2NO3_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SiOX = IndForcing('SiOH', 'Forcing/constantMLD/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/MLDdriven/SST_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/MLDdriven/PAR_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.type = 'flowthrough'

        elif forcingtype == 'batch':
            self.MLD = IndForcing('MLD', 'Forcing/MLDdriven/MLD2015_R1.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.NOX = IndForcing('NO2NO3', 'Forcing/constantMLD/NO2NO3_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SiOX = IndForcing('SiOH', 'Forcing/constantMLD/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/MLDdriven/SST_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/MLDdriven/PAR_R1.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.type = 'batch'

        #WOA STATIONS EMPOWER:
        elif forcingtype == 'BIOTRANS':
            self.MLD = IndForcing('mld', 'X', k=5, s=100, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=47, Lon=-20, RBB=2.5)
            self.NOX = IndForcing('n0x', 'X', k=5, s=None, kind="PWPoly", forctype=forcingtype,
                                  WOA=True, Lat=47, Lon=-20, RBB=2.5)
            self.SiOX = IndForcing('SiOH', 'Forcing/SiOH_filler/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('sst', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=47, Lon=-20, RBB=2.5)
            self.PAR = IndForcing('par', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=47, Lon=-20, RBB=2.5)
            self.verif = VerifData(Lat=47, Lon=-20, RBB=2.5)
            self.type = 'MLD'

        elif forcingtype == 'PAPA':
            self.MLD = IndForcing('mld', 'X', k=5, s=100, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=50.1, Lon=-144.9, RBB=2.5)
            self.NOX = IndForcing('n0x', 'X', k=5, s=None, kind="PWPoly", forctype=forcingtype,
                                  WOA=True, Lat=50.1, Lon=-144.9, RBB=2.5)
            self.SiOX = IndForcing('SiOH', 'Forcing/SiOH_filler/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('sst', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=50.1, Lon=-144.9, RBB=2.5)
            self.PAR = IndForcing('par', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=50.1, Lon=-144.9, RBB=2.5)
            self.verif = VerifData(Lat=50.1, Lon=-144.9, RBB=2.5)
            self.type = 'MLD'

        elif forcingtype == 'CARIACO':
            self.MLD = IndForcing('mld', 'X', k=5, s=100, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=10.5, Lon=-64.6, RBB=2.5)
            self.NOX = IndForcing('n0x', 'X', k=5, s=None, kind="PWPoly", forctype=forcingtype,
                                  WOA=True, Lat=10.5, Lon=-64.6, RBB=2.5)
            self.SiOX = IndForcing('SiOH', 'Forcing/SiOH_filler/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('sst', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=10.5, Lon=-64.6, RBB=2.5)
            self.PAR = IndForcing('par', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=10.5, Lon=-64.6, RBB=2.5)
            self.type = 'MLD'

        elif forcingtype == 'BATS':

            self.MLD = IndForcing('mld', 'X', k=5, s=100, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=31, Lon=-64, RBB=2.5)
            self.NOX = IndForcing('n0x', 'X', k=5, s=None, kind="PWPoly", forctype=forcingtype,
                                  WOA=True, Lat=31, Lon=-64, RBB=2.5)
            self.SiOX = IndForcing('SiOH', 'Forcing/SiOH_filler/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('sst', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=31, Lon=-64, RBB=2.5)
            self.PAR = IndForcing('par', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=31, Lon=-64, RBB=2.5)
            self.verif = VerifData(Lat=31, Lon=-64, RBB=2.5)
            self.type = 'MLD'

        elif forcingtype == 'WOA2018':

            self.MLD = IndForcing('mld', 'X', k=5, s=100, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=lat, Lon=lon, RBB=rbb)
            self.NOX = IndForcing('n0x', 'X', k=5, s=None, kind="PWPoly", forctype=forcingtype,
                                  WOA=True, Lat=lat, Lon=lon, RBB=rbb)
            self.SiOX = IndForcing('SiOH', 'Forcing/SiOH_filler/SiOH_R1.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('sst', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=lat, Lon=lon, RBB=rbb)
            self.PAR = IndForcing('par', 'X', k=5, s=5, kind="spline", forctype=forcingtype,
                                  WOA=True, Lat=lat, Lon=lon, RBB=rbb)
            self.verif = VerifData( Lat=lat, Lon=lon, RBB=rbb)
            self.type = 'MLD'

        elif forcingtype == 'EMPOWER':
            self.MLD = IndForcing('MLD', 'Forcing/EMPOWER/EMPOWER-forcing-export.csv', k=5, s=100, kind="spline", forctype=forcingtype)
            self.X21 = IndForcing('depth', 'Forcing/X21Iso/X21Iso_r1.csv', k=5, s=100, kind="spline",
                                  forctype=forcingtype)
            self.NOX = IndForcing('N0', 'Forcing/EMPOWER/EMPOWER-forcing-export.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.SST = IndForcing('SST', 'Forcing/EMPOWER/EMPOWER-forcing-export.csv', k=5, s=None, kind="PWPoly", forctype=forcingtype)
            self.PAR = IndForcing('PAR', 'Forcing/EMPOWER/EMPOWER-forcing-export.csv', k=5, s=None, kind="spline", forctype=forcingtype)
            self.verif = VerifData(Lat=47, Lon=-20, RBB=2.5)
            self.type = 'MLD'

        else:
            raise('wrong forcingtype passed to Forcing class')

class Plankton:
    """
    initializes the Plankton (both Zoo- & Phyto-) community according to number of types prescribed by params
    """

    def __init__(self, modelparams, zptype):
        if zptype == 'Zooplankton':
            self.zoonum = modelparams['zoo_num'].value
            self.stdparams = sliceoffparams(modelparams, 'zt')
            self.pars = [ZooType(self.stdparams, sliceparams(modelparams, 'zt' + str(i + 1)))
                        for i in range(self.zoonum)]
        elif zptype == 'Phytoplankton':
            self.phytonum = modelparams['pfun_num'].value
            self.stdparams = sliceoffparams(modelparams, 'pt')
            self.pars = [PhytoType(self.stdparams, sliceparams(modelparams, 'pt' + str(i + 1))) for i in
                         range(self.phytonum)]

        else:
            raise('wrong functional type passed to Plankton class')

    def init(self, ):
        return self.pars

class PhytoType:
    """sets up individual functional types and defines functions
    such as nutrient uptake and growth"""

    def __init__(self, stdpars, slicedpars):
        self.kw = checkreplaceparam(stdpars, slicedpars, 'kw')
        self.OptI = checkreplaceparam(stdpars, slicedpars, 'OptI')

        self.U_N = checkreplaceparam(stdpars, slicedpars, 'U_N')
        self.U_Si = checkreplaceparam(stdpars, slicedpars, 'U_Si')

        self.v = checkreplaceparam(stdpars, slicedpars, 'v')

        self.muP = checkreplaceparam(stdpars, slicedpars, 'muP')
        self.moP = checkreplaceparam(stdpars, slicedpars, 'moP')

        self.ratioSi = checkreplaceparam(stdpars, slicedpars, 'ratioSi')

        self.pfn = stdpars['pfun_num']
        self.zn = stdpars['zoo_num']

        # print('PHYTO',self.pfn, self.zn)

        self.kc = checkreplaceparam(stdpars, slicedpars, 'kc')
        self.alpha = checkreplaceparam(stdpars, slicedpars, 'alpha')
        self.VpMax = checkreplaceparam(stdpars, slicedpars, 'VpMax')

        self.zoolist = ['Z'+str(j+1) for j in range(self.zn)]
        # print(self.zoolist)
        self.grazepref = [checkreplaceparam(stdpars, slicedpars, string) for string in self.zoolist]

        # print('grazepref', self.grazepref)


    def n_uptake(self, Nitrate):
        N_Uptake = Nitrate / (Nitrate + self.U_N)  # Michaelis Menten - uptake of Nitrate
        return N_Uptake

    def si_uptake(self, Silicate):
        if self.U_Si == 0:
            # non-diatoms
            return 0
        else:
            # diatoms
            Si_Uptake = Silicate / (Silicate + self.U_Si)  # Michaelis Menten - uptake of Nitrate
            return Si_Uptake

    def lightharvesting(self, intMLD, intPAR):
        # light harvesting according to Steele PI curve, and a constant extinction coefficient
        self.kw = 0.1
        lighthrv =  1. / (self.kw * intMLD) * \
               (-np.exp(1. - intPAR / self.OptI) - (-np.exp((1. - (intPAR * np.exp(-self.kw * intMLD)) / self.OptI))))
        return lighthrv

    def smithpi(self, intMLD, intPAR, P):
        # lightharvesting according to the smith equation (from EMPOWER publication)
        # still need to add kPAR calculation based on kw and kchla
        kPAR = self.kw + self.kc * sum(P)
        x_0 = self.alpha * intPAR * np.exp(- kPAR * 0)
        x_H = self.alpha * intPAR * np.exp(- kPAR * intMLD)
        VpH = (self.VpMax / (kPAR * intMLD)) * \
              np.log(
                  (x_0 + np.sqrt(self.VpMax ** 2 + x_0 ** 2)) /
                  (x_H + np.sqrt(self.VpMax ** 2 + x_H ** 2))
              )
        return VpH

    def tempdepgrowth(self, int_SST):
        tdp = np.exp(0.063 * int_SST) ##############################################################################
        return tdp

    def gains(self, nuptake, siuptake, lighthrv, tempdepgro, P):
        if self.U_Si == 0:
            # non-diatoms
            Gains = (self.muP * nuptake * lighthrv * tempdepgro) * P
        else:
            # diatoms
            Gains = (self.muP * min(nuptake, siuptake) * lighthrv * tempdepgro) * P
        return Gains

    def sinking(self, intMLD, P):
        Sink = self.v / intMLD * P  # Phytoplankton sinking as a function of MLD and sinking rate
        return Sink

    def silicatedrawdown(self, Gainz):
        if self.U_Si == 0:
            # non-diatoms
            return 0
        else:
            # diatoms
            Si_Drawdown = self.ratioSi * Gainz
            return Si_Drawdown

    def mortality(self, P):
        Mortal = self.moP * P
        return Mortal

    def zoograzing(self, Gj, Pi, Z):
        # take the general grazing term from each zooplankton, multiply by phyto fraction and sum
        Grazing = [Gj[j] * (self.grazepref[j] * Pi ** 2) * Z[j] for j in range(self.zn)]
        GrazingPerZ = sum(Grazing)
        return GrazingPerZ



class ZooType:
    def __init__(self, stdpars, slicedpars):
        # zooplankton
        self.moZ = checkreplaceparam(stdpars, slicedpars, 'moZ')
        self.muZ = checkreplaceparam(stdpars, slicedpars, 'muZ')
        # grazing params
        self.Kp = checkreplaceparam(stdpars, slicedpars, 'Kp')
        self.deltaZ = checkreplaceparam(stdpars, slicedpars, 'deltaZ')

        self.muIntGraze = checkreplaceparam(stdpars, slicedpars, 'muIntGraze')
        self.kIntGraze = checkreplaceparam(stdpars, slicedpars, 'kIntGraze')

        self.pred = checkreplaceparam(stdpars, slicedpars, 'pred')
        self.deltaLambda = checkreplaceparam(stdpars, slicedpars, 'deltaLambda')

        self.pfn = stdpars['pfun_num']
        self.zn = stdpars['zoo_num']
        # print('ZOO', self.pfn,self.zn)

        self.beta = 2 # for Vallina KTW Grazing, beta = 1 : Holling Type II, beta = 2 : Holling Type III
        self.ksat = self.Kp

        self.phylist = ['P' + str(i + 1) for i in range(self.pfn)]
        self.zoointlistfeed = ['Zint_feed' + str(j + 1) for j in range(self.zn)] # list of feeding
        self.zoointlistgrazed = ['Zint_grazed' + str(j + 1) for j in range(self.zn)]
        # print(self.phylist, self.zoointlistfeed, self.zoointlistgrazed)

        self.feedpref = [checkreplaceparam(stdpars, slicedpars, string) for string in self.phylist]
        self.interfeedpref = [checkreplaceparam(stdpars, slicedpars, string) for string in self.zoointlistfeed]
        self.intergrazedpref = [checkreplaceparam(stdpars, slicedpars, string) for string in self.zoointlistgrazed]

        # print('feedpref',self.feedpref, 'intfeedpref', self.interfeedpref, 'intrgrpref', self.intergrazedpref, 'muZ', self.muZ)

    def zoomortality(self, Z):
        # linear mortality of zooplankton, quadratic loss via higher order pred (below)
        total_moZ = self.moZ * Z
        return total_moZ

    def zoofeeding(self, P, Z, func='anderson'):
        if func == 'anderson':
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP + FrhoZ
            GrazingProb = self.muZ / (self.ksat ** 2 + Frho)
            return GrazingProb

        elif func == 'fasham': #put the holling type num here
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP + FrhoZ  # active switching coefficient
            FpP = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])
            FpZ = sum([self.interfeedpref[j] * Z[j] for j in range(self.zn)])
            Fp = FpP + FpZ  # total food available
            GrazingProb = self.muZ * (1 / (self.ksat * Fp + Frho))
            return GrazingProb

        elif func == 'vallina':
            FrhoP = sum([self.feedpref[i] * P[i] ** 2 for i in range(self.pfn)])
            FrhoZ = sum([self.interfeedpref[j] * Z[j] ** 2 for j in range(self.zn)])
            Frho = FrhoP + FrhoZ  # active switching coefficient
            FpP = sum([self.feedpref[i] * P[i] for i in range(self.pfn)])
            FpZ = sum([self.interfeedpref[j] * Z[j] for j in range(self.zn)])
            Fp = FpP + FpZ  # total food available
            GrazingProb = self.muZ * (1 / Frho) * ((Fp ** self.beta) / (self.ksat ** self.beta + Fp ** self.beta))
            return GrazingProb

        else:
            print('no grazing formulation given, wrong func key')

    def fullgrazing(self, GjJ, P, Z, Zi):
        # phytouptake + zooplankton per zooplankton for each phyto
        IprobP = [GjJ * (self.feedpref[i] * P[i] ** 2) for i in range(self.pfn)]  # grazeprob per each PFT
        IprobZ = [GjJ * (self.interfeedpref[j] * Z[j] ** 2) for j in range(self.zn)]
        Iprob = IprobP + IprobZ
        #perhaps here = multiply by Z before summing Iprobs!
        Itots = sum(Iprob)
        Itot = Itots * Zi
        return Itot

    def assimgrazing(self, Itoti):
        AssimGrazing = self.deltaZ * Itoti
        return AssimGrazing

    def unassimilatedgrazing(self, Itoti):
        UnAsGraze = (1. - self.deltaZ) * Itoti
        return UnAsGraze

    def interzoograze(self, Gj, Z, Zi):
        # returns total biomass consumed by interzoopred per zoo type (being grazed upon)
        IntGraze = [Gj[j] * (self.intergrazedpref[j] * Zi ** 2) * Z[j] for j in range(self.zn)]
        TotIntGraze = sum(IntGraze)
        return TotIntGraze

    def higherorderpred(self, Zi):
        total_moZ = self.pred * Zi ** 2  # sum(Z)
        return total_moZ / self.zn
