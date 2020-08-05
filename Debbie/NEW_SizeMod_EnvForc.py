import os
import numpy as np
import scipy.interpolate as intrp

# Import forcing data
#os.chdir('/Users/szewing/Desktop/PhD_work/SizeMod/Forcing/Data')
# os.chdir('/people/home/sto/SizeMod/Forcing')


class SizeMod_EnvForc:
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
        self.Zmix = Zmix                                      # mixing depth [meter]
        self.ThermoZ = ThermoZ                                # thermocline depth [meter]
        
        self.data = self.data_extract()
        print(self.fx_type, self.data, type(self.data))
        self.interpolated = self.intrp_forcing(self.data)     # storing interpolated forcing
        self.deriv = self.interpolated.derivative()                # storing interpolated derivative



    def nSSI2PAR(self, netSSI):
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
        μmol2mol = 1e6
        coef = extr * conv * sec2day / μmol2mol  # 0.1776816
        return netSSI * coef

    def data_extract(self):
        if self.fx_type == 'LWST':
            data_name = 'LWST_13p.csv'
            LWST_data = np.genfromtxt(data_name, delimiter=',')
            LWST_data = LWST_data[:, 5]  # extracting temperature data of 40N
            return LWST_data
        elif self.fx_type == 'PAR':
            data_name = 'nSSI_13p.csv'
            PAR_data = np.genfromtxt(data_name, delimiter=',')
            PAR_data = PAR_data[:, 5]   # extracting irradiance data of 40N
            PAR_data = self.nSSI2PAR(PAR_data)
            return PAR_data


    def return_intrp(self, timestep):
        """
        Function for returning interpolated value at time step
        :param timestep:
        """
        NewT = np.mod(timestep, 365.)
        return self.interpolated(NewT)


    def return_deriv(self, timestep):
        """
        Function for returning derivative from interpolated forcing
        :param timestep:
        """
        NewT = np.mod(timestep, 365.)
        return self.deriv(NewT)


    def intrp_forcing(self, data):
        """
        This function perform interpolation for all environmental forcing (temperature, irradiance and MLD)
        :param data: should be in np array format. Only for temperature and par, for mixing no data is required
        :return: an interpolate object spl
        """
        #NewT = np.mod(timestep, 365.)  # to convert 'time' into day of the year after the first year
        x = np.arange(365)

        if self.fx_type == 'LWST':
            spl = intrp.CubicSpline(np.linspace(0,365, len(data)+1), np.concatenate([data,data[-1]], axis=None))
            #lwst = spl(NewT)
            return spl
        elif self.fx_type == 'PAR':
            spl = intrp.CubicSpline(np.linspace(0,365, len(data)+1), np.concatenate([data,data[-1]], axis=None))
            #par = spl(NewT)
            return spl
        elif self.fx_type == 'MLD' or 'dMdt':
            if self.freq == 'constant':
                spl = intrp.CubicSpline(x, np.full(len(x), self.Zmix))
                return spl
            elif self.freq == 'medium':
                mld = (-((self.Zmix - self.ThermoZ) / 2 + self.ThermoZ) + (self.Zmix - self.ThermoZ) / 2 *
                       np.cos(x / 14.525)) * -1
                spl = intrp.CubicSpline(x, mld)
                return spl
            elif self.freq == 'high':
                mld = (-((self.Zmix - self.ThermoZ) / 2 + self.ThermoZ) + (self.Zmix - self.ThermoZ) / 2 *
                       np.cos(x / 4.825)) * -1
                spl = intrp.CubicSpline(x, mld)
                return spl

