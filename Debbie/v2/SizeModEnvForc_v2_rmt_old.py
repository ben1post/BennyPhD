import os
import numpy as np
import scipy.interpolate as intrp

# Import forcing data
#os.chdir('/Users/szewing/Desktop/PhD_work/SizeMod/Forcing/Data')
#os.chdir('/people/home/sto/SizeMod/Forcing')


class SizeModEnvForc:
    """
    Environmental forcing for the size-based model
    Method to interpolate from monthly to daily environmental forcing data, for:
    (i)     Mixed layer depth (MLD);
    (ii)    Photosynthetically active radiation (PAR)
    (iii)   Lake water surface temperature (LWST)
    (iv)    dMLDdt
    ----------
    The interpolation is determined by cubic spline approximation,
    this method interpolates data with a piecewise cubic polynomial which is twice continuously differentiable
    (ref: https://docs.scipy.org/doc/scipy/reference/interpolate.html)
    ----------
    https://math.stackexchange.com/questions/100655/cosine-esque-function-with-flat-peaks-and-valleys
    The function produces interpolation for the theoretical sinusoidal mixing depth
    """

    def __init__(self, fx_type, freq=None, Zmix=80, ThermoZ=2.5):
        """
        :param fx_type: LWST, PAR, MLD
        :param freq: constant, medium, high
        :param Zmix: Mixing depth (meter)
        :param ThermoZ: Thermocline depth (meter)
        """
        self.fx_type = fx_type
        self.freq = freq
        self.Zmix = Zmix        # mixing depth [meter]
        self.ThermoZ = ThermoZ  # thermocline depth [meter]

    @staticmethod
    def nSSI2PAR(netSSI):
        """
        This function is to convert nSSI (W m-2) into PAR (Ein m-2 d-1/mol m-2 d-1)
        PAR is a term used to describle radiation in wavelengths usedful for photosynthesis [1]
        PAR is generally accepted to be wavelengths between 400 and 700 nm [1]
        Procedure:
        1. extr = extraction factor (only about 45% of the nSSI is in the 400-700nm range)
        2. conv = W m-2 converts into μmole m-2 s-1 = 4.57 [1,2]
        3. sec2day = converting from s-1 to d-1
        4. μmol2mol = converting from μmol to mol (1 mole = 1 Ein)

        [1] Sager and Farlane, 1997
        Sager, J.C., McFarlane, J.C. 1997. Radiation. In: Plant Growth Chamber Handbook. (Langhans, R.W. and
        Tibbitts, T.W., Eds.) Iowa State University: North Central Regional Research Publication No. 340, Iowa
        Agriculture and Home Economics Experiment Station Special Report no. 99, pp. 1-29.
        [2] Thimijan and Heins, 1983
        Thimijan, R.W., and R.D. Heins. 1983. Photometric, radiometric and quantum light units of measure:
        A review of procedures for interconversion. HortScience 18(6): 818-822
        """
        extr = 0.45
        conv = 4.57
        sec2day = 86400
        umol2mol = 1e6
        coef = extr * conv * sec2day / umol2mol  # 0.1776816
        return netSSI * coef


    def intrp_forcing(self, tstep):
        """
        This function perform interpolation for all environmental forcing (temperature, irradiance and MLD)
        :param Ynum: number of years
        :param tstep: time step
        :return: an interpolate object spl
        """

        NewT = np.mod(tstep, 365.)  # to convert 'time' into day of the year after the first year
        x = np.arange(365)
        # c1 = 2
        # c3 = -1
        dayspermonth = [0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        dayofyear = np.cumsum(dayspermonth)

        if self.fx_type == 'LWST':
            data_name = 'LWST_13p_syn.csv'  # set up data name for import
            LWST_data = np.genfromtxt(data_name, delimiter=',')
            spl_t = LWST_data[:, 0]
            #spl_t = dayofyear
            LWST_data = LWST_data[:, 5]  # extracting data from 40N
            data = LWST_data
            spl = intrp.UnivariateSpline(spl_t, data, k=4, s=1)
            return spl(NewT)
        elif self.fx_type == 'PAR':
            data_name = 'nSSI_13p_syn.csv'
            PAR_data = np.genfromtxt(data_name, delimiter=',')
            spl_t = PAR_data[:, 0]
            #spl_t = dayofyear
            PAR_data = PAR_data[:, 5]
            PAR_data = self.nSSI2PAR(PAR_data)
            data = PAR_data
            spl = intrp.UnivariateSpline(spl_t, data, k=4, s=1)
            return spl(NewT)
        elif self.fx_type == 'MLD':     # or 'dMdt':
            if self.freq == 'constant':
                data = np.full(len(x), self.Zmix)
                spl = intrp.UnivariateSpline(x, data, k=4, s=1)
                out = spl(NewT)
            elif self.freq == 'medium':
                out = (-((self.Zmix - self.ThermoZ) / 2 + self.ThermoZ) + (self.Zmix - self.ThermoZ) / 2 *
                             np.cos(NewT / 14.525)) * -1
            elif self.freq == 'high':
                out = (-((self.Zmix - self.ThermoZ) / 2 + self.ThermoZ) + (self.Zmix - self.ThermoZ) / 2 *
                             np.cos(NewT / 4.825)) * -1
            return out
        elif self.fx_type == 'dMdt':
            if self.freq == 'constant':
                data = np.full(len(x), self.Zmix)
                spl = intrp.UnivariateSpline(x, data, k=4, s=1)
                out = spl.derivative()(NewT)
            elif self.freq == 'medium':
                out = 1 * (-self.ThermoZ + self.Zmix) * np.sin(NewT / 14.525) / (2 * 14.525)
                # derivative function obtained from sympy
            elif self.freq == 'high':
                out = 1 * (-self.ThermoZ + self.Zmix) * np.sin(NewT / 4.825) / (2 * 4.825)
            return out
