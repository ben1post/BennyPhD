import pandas
import numpy as np
import scipy.interpolate as intrp
#from phydra.aux import sliceparams, sliceoffparams, checkreplaceparam

from scipy.io import netcdf
import os

from runs.modelrun_fullTS import parameters
from phydra.core import ModelSetup

ms = ModelSetup(parameters, physics='Box', forcing='fullTS', time=None)

n, p, z, d = ms.classes

physx = ms.physics

print(physx.forcing.NOX)


