#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas
import numpy as np
import scipy.interpolate as intrp

from scipy.io import netcdf
import os

# TODO:
#  - add missing nutrients
#  - add HPLC functype data

class CARIACOdata:
    """
    initializes and reads new forcing & verification data of full time series & supplies this to model
    either aggregated per regime or as full time series
    """
    def __init__(self, forcvar, data='niskin', time='regime1', k=3, s=None, kind="spline", boxordep='box', forctype='aggTS', pad=False, extendstart=False):
        """# parameters for interpolation"""
        self.k = k
        self.s = s
        self.kind = kind
        self.forcingtype = forctype
        self.forcvar = forcvar
        self.pad = pad
        self.extendstart = extendstart
        self.time = time

        if data == 'niskin':
            if boxordep == 'box':
                varname = forcvar + '_Box'
            elif boxordep == 'depth':
                varname = forcvar + '_AtDepth'
            self.varname = varname
        else:
            varname = forcvar
            self.varname = varname

        self.forcingfile = self.readCariaco(forcvar, data=data, time=time, boxordep=boxordep, forctype=forctype)
        self.rawforcing = self.readCariaco(forcvar, data=data, time=time, boxordep=boxordep, forctype=forctype, get='full')

        if forctype == 'aggTS':
            self.interpolated = self.dailyinterp(self.forcingfile, self.kind, self.k, self.s)
            if kind == "spline":
                self.derivative = self.interpolated.derivative()
                self.derivative2 = self.interpolated.derivative(n=2)

        elif forctype == 'fullTS':
            self.interpolated_df = self.fullinterp(varname=forcvar, data=data, boxordep=boxordep)
            self.interpolated = self.dailyinterp(self.forcingfile, self.kind, self.k, self.s)
            if kind == "spline":
                self.derivative = self.interpolated.derivative()

        print(forcvar + ' forcing created')

    def returnModelOut_nospinup(self, outarray):
        if self.forcingtype == 'fullTS':
            if self.pad:
                if self.extendstart:
                    #8767 vs 7737
                    return outarray[8767-7737:8767]
        print('No need to remove spinup phase')
        return outarray

    def returnTimeDays_nospinup(self, timedays):
        if self.forcingtype == 'fullTS':
            if self.pad:
                if self.extendstart:
                    #8767 vs 7737
                    return timedays[:7737]
        print('No need to remove spinup phase')
        return timedays

    def returnModelOut_Regime(self, outarray, timedays, regime, spinup=False):
        if not spinup:
            if regime == 1:
                b4_reg1_df = self.interpolated_df.loc[:'1996-01-01']
                reg1_df = self.interpolated_df.loc['1996-01-01':'2000-10-30']
                a = len(b4_reg1_df)
                b = len(reg1_df)
                out = outarray[a:a+b]
                time = timedays[:b]
                return [out,time]
            elif regime == 2:
                b4_reg2_df = self.interpolated_df.loc[:'2006-01-01']
                reg2_df = self.interpolated_df.loc['2006-01-01':'2010-12-31']
                a = len(b4_reg2_df)
                b = len(reg2_df)
                out = outarray[a:a+b]
                time = timedays[:b]
                return [out, time]


    def readCariaco(self, varname, data, time, boxordep, forctype, get=None):
        """read data and return either aggregated regimes or full time series"""
        if data == 'niskin':
            df_all = pandas.read_csv('Data/NewestData/BoxVSatDepth_02.csv')
            if boxordep == 'box':
                varname = varname + '_Box'
            elif boxordep == 'depth':
                varname = varname + '_AtDepth'
                df_all[varname] = np.nanmean(df_all[varname])
        elif data == 'SeaWiFS':
            df_all = pandas.read_csv('Data/NewestData/PARXSeaWiFS_03.csv')
        elif data == 'x25.8':
            df_all = pandas.read_csv('Data/NewestData/x258_02.csv')

        df_all.date = pandas.to_datetime(df_all.date)
        df_val = df_all[['date', 'month', 'yday', varname]]
        df_series = df_val.set_index(df_val['date'])

        if self.extendstart:
            new_row = pandas.DataFrame({'date':pandas.to_datetime(['1993-01-01']), 'month': 1, 'yday': 1, varname: np.NaN})
            new_row = new_row.set_index(new_row['date'])

            df_X = pandas.concat([new_row, df_series])

            df_X_M = df_X.resample('M').mean()
            df_series = df_X_M


            if self.forcvar == 'x258depth' and self.pad == False:

                data_pd = df_series
                data_pd.fillna(data_pd.groupby(data_pd.index.month).transform('mean'), inplace=True)

                df_series_cut = data_pd.loc['1993-01-01':'1995-11-08']
                df_series_rest = df_series.loc['1995-11-08':]
                df_series = pandas.concat([df_series_cut, df_series_rest])

        if self.pad:
            data_pd = df_series
            data_pd.fillna(data_pd.groupby(data_pd.index.month).transform('mean'), inplace=True)
            df_series = data_pd



        if self.forcingtype == 'aggTS':
            if self.time == 'regime1':
                df_series_cut = df_series.loc['1996-01-01':'2000-10-30']
                df_val = df_series_cut
            elif self.time == 'regime2':
                df_series_cut = df_series.loc['2006-06-30':'2010-12-31']
                df_val = df_series_cut
            if get == 'full':
                return df_val
            df_monthly_mean = df_val.groupby('month').mean()
            forcing_oneyear = list(df_monthly_mean[varname])
            forcing_list = forcing_oneyear * 3
            return forcing_list

        elif self.forcingtype == 'fullTS':
            df_val = df_series
            return df_val

    def fullinterp(self, data, varname, boxordep):
        """test"""
        if data == 'niskin':
            if boxordep == 'box':
                varname = varname + '_Box'
            elif boxordep == 'depth':
                varname = varname + '_AtDepth'

        df = self.forcingfile.copy()

        df_days = df.resample('D').mean()
        df_days[varname] = df_days[varname].interpolate()
        return df_days

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
        if self.forcingtype == 'aggTS':
            dayspermonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            dpm = dayspermonth * 3
            dpm_cumsum = np.cumsum(dpm) - np.array(dpm) / 2
            if kind == 'spline':
                outintp = intrp.UnivariateSpline(dpm_cumsum, file, k=k, s=s)
                return outintp
            elif kind == 'PWPoly':
                outintp = intrp.PchipInterpolator(dpm_cumsum, file)
                return outintp
        elif self.forcingtype == 'fullTS':
            numrows = np.arange(0,self.interpolated_df.shape[0],1)
            #print('fix me!',numrows, self.interpolated_df[self.varname].values)
            outintp = intrp.UnivariateSpline(numrows, np.array(self.interpolated_df[self.varname].values), k=self.k, s=self.s)
            return outintp
        else:
            raise ('Wrong interpolation type passed to dailyinterp function of IndForcing class')


    def return_interpvalattime(self, t):
        if self.forcingtype == 'fullTS':
            return self.interpolated(t)
        else:
            """
            Method to return interpolated value of forcing.

            converts time in days to time in months
            """
            newt = np.mod(t, 365.) + 365  # *12./365.
            return self.interpolated(newt)


    def return_derivattime(self, t):
        """
        Method to return derivative (slope) of interpolated value of forcing.

        converts time in days to time in months, and the resulting derivative from per month to per day
        """
        if self.forcingtype == 'fullTS':
            return self.derivative(t)

        newt = np.mod(t, 365.) + 365  # * 12. / 365.

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
    def __init__(self, forcingtype, time, pad, extendstart):
        if forcingtype == 'aggTS':
            self.X258 = CARIACOdata('x258depth', data='x25.8', time=time, k=3, s=2000, kind="spline", forctype=forcingtype, pad=pad)
            self.NOX = CARIACOdata('NO3_NO2_USF', data='niskin', time=time, boxordep='depth', k=50, s=None, kind="PWPoly", forctype=forcingtype, pad=pad)
            self.SiOH = CARIACOdata('SiO4_USF', data='niskin', time=time, boxordep='depth', k=50, s=None, kind="PWPoly", forctype=forcingtype, pad=pad)
            self.SST = CARIACOdata('Temperature', data='niskin', time=time, k=2, s=5, kind="spline", forctype=forcingtype, pad=pad)
            self.PAR = CARIACOdata('PAR', data='SeaWiFS', time=time, k=5, s=None, kind="spline", forctype=forcingtype, pad=pad)
            self.verif = CARIACOVerifData(time=time, forctype=forcingtype, pad=pad)
            self.type = 'Box'

        elif forcingtype == 'fullTS':
            self.X258 = CARIACOdata('x258depth', data='x25.8', time=time, k=3, s=2000000, kind="spline", forctype=forcingtype, pad=True, extendstart=extendstart)
            self.NOX = CARIACOdata('NO3_NO2_USF', data='niskin', time=time, boxordep='depth', k=2, s=100, kind="PWPoly", forctype=forcingtype, pad=pad, extendstart=extendstart)
            self.SiOH = CARIACOdata('SiO4_USF', data='niskin', time=time, boxordep='depth', k=2, s=100, kind="PWPoly", forctype=forcingtype, pad=pad, extendstart=extendstart)
            self.SST = CARIACOdata('Temperature', data='niskin', time=time, k=2, s=5, kind="spline", forctype=forcingtype, pad=pad, extendstart=extendstart)
            self.PAR = CARIACOdata('PAR', data='SeaWiFS', time=time, k=5, s=None, kind="spline", forctype=forcingtype, pad=pad, extendstart=extendstart)
            self.verif = CARIACOVerifData(time=time, forctype=forcingtype, pad=pad)
            self.type = 'Box'

        else:
            raise('wrong forcingtype passed to Forcing class')


class CARIACOVerifData:
    """initializes and reads verification data from CARIACO time series"""
    def __init__(self, time, forctype, pad):
        self.pad = False
        self.time = time
        self.forcingtype = forctype
        self.HPLC = self.readHPLC('Tchla')
        self.FluorChla = self.readniskin('Chlorophyll', boxordep='box')
        self.NO2NO3 = self.readniskin('NO3_NO2_USF', boxordep='box')
        self.SiOH = self.readniskin('SiO4_USF', boxordep='box')
        self.PN = self.readniskin('PON_ug_kg', boxordep='box')
        self.Zoo = self.readZoo('BIOMASS_200')
        # self.pad = pad
        self.padHPLC = self.readHPLC('Tchla')
        self.fullpadHPLC = self.fullinterp(self.padHPLC, 'Tchla')

        self.padFluorChla = self.readniskin('Chlorophyll', boxordep='box')
        self.fullpadFluorChla = self.fullinterp(self.padFluorChla, 'Chlorophyll', boxordep='box', data='niskin')

        self.padNO2NO3 = self.readniskin('NO3_NO2_USF', boxordep='box')
        self.fullpadNO2NO3 = self.fullinterp(self.padNO2NO3, 'NO3_NO2_USF', boxordep='box', data='niskin')

        self.padSiOH = self.readniskin('SiO4_USF', boxordep='box')
        self.fullpadSiOH = self.fullinterp(self.padSiOH, 'SiO4_USF', boxordep='box', data='niskin')

        self.padPN = self.readniskin('PON_ug_kg', boxordep='box')
        self.fullpadPN = self.fullinterp(self.padPN, 'PON_ug_kg', boxordep='box', data='niskin')

        self.padZoo = self.readZoo('BIOMASS_200')
        self.fullpadZoo = self.fullinterp(self.padZoo, 'BIOMASS_200')

        print('Verification data created')

    def returnMeanVerifPerMonth(self, dataset, varname):
        if self.forcingtype == 'fullTS':
            return dataset
        df_monthly_mean = dataset.groupby('month').mean()
        verif_oneyear = list(df_monthly_mean[varname])
        return verif_oneyear

    def readHPLC(self, varname):
        df_all = pandas.read_csv('Data/NewestData/HPLCPinckneyTotAndSpec_03.csv')
        df_all.date = pandas.to_datetime(df_all.date)
        try:
            df_val = df_all[['date', 'month', 'yday', varname]]
        except:
            print(df_all.columns)
            raise Exception('Variable {} is not found in HPLC dataframe'.format(varname))

        df_series = df_val.set_index(df_val['date'])

        if self.pad:
            data_pd = df_series
            data_pd.fillna(data_pd.groupby(data_pd.month).transform('mean'), inplace=True)
            df_series = data_pd

        if self.forcingtype == 'aggTS':
            if self.time == 'regime1':
                df_series_cut = df_series.loc['1996-01-01':'2000-10-30']
                df_val = df_series_cut
            elif self.time == 'regime2':
                df_series_cut = df_series.loc['2006-06-30':'2010-12-31']
                df_val = df_series_cut
        if self.forcingtype == 'fullTS':
            df_val = df_series

        return df_val


    def readniskin(self, varname, boxordep):
        df_all = pandas.read_csv('Data/NewestData/BoxVSatDepth_02.csv')
        df_all.date = pandas.to_datetime(df_all.date)

        if boxordep == 'box':
            varname = varname + '_Box'
        elif boxordep == 'depth':
            varname = varname + '_AtDepth'
        try:
            df_val = df_all[['date', 'month', 'yday', varname]]
        except:
            print(df_all.columns)
            raise Exception('Variable {} is not found in niskin df'.format(varname))

        df_series = df_val.set_index(df_val['date'])

        if self.pad:
            data_pd = df_series
            data_pd.fillna(data_pd.groupby(data_pd.month).transform('mean'), inplace=True)
            df_series = data_pd

        if self.forcingtype == 'aggTS':
            if self.time == 'regime1':
                df_series_cut = df_series.loc['1996-01-01':'2000-10-30']
                df_val = df_series_cut
            elif self.time == 'regime2':
                df_series_cut = df_series.loc['2006-06-30':'2010-12-31']
                df_val = df_series_cut
        if self.forcingtype == 'fullTS':
            df_val = df_series
        return df_val

    def readZoo(self, varname):
        df_all = pandas.read_csv('Data/NewestData/ZooplanktonData_05.csv')
        df_all.date = pandas.to_datetime(df_all.date)
        try:
            df_val = df_all[['date', 'month', 'yday', 'BIOMASS_200', 'BIOMASS_500']]
        except:
            print(df_all.columns)
            raise Exception('Variable {} is not found in ZOO dataframe'.format(varname))
        df_series = df_val.set_index(df_val['date'])

        if self.pad:
            data_pd = df_series
            data_pd.fillna(data_pd.groupby(data_pd.month).transform('mean'), inplace=True)
            df_series = data_pd

        if self.forcingtype == 'aggTS':
            if self.time == 'regime1':
                #print(df_series.head(10))
                df_series_cut = df_series.loc['1996-01-01':'2000-10-30']
                df_val = df_series_cut
            elif self.time == 'regime2':
                df_series_cut = df_series.loc['2006-06-30':'2010-12-31']
                df_val = df_series_cut

        if self.forcingtype == 'fullTS':
            df_val = df_series

        return df_val

    def fullinterp(self, df, varname, boxordep=None, data=None):
        """test"""
        if data == 'niskin':
            if boxordep == 'box':
                varname = varname + '_Box'
            elif boxordep == 'depth':
                varname = varname + '_AtDepth'

        df_days = df.resample('D').mean()

        #df_days[varname] = df_days[varname].interpolate()
        return df_days
