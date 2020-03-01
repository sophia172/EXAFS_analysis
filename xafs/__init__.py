#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,UnivariateSpline,LSQUnivariateSpline
from scipy.optimize import curve_fit
import glob
from scipy import stats
from scipy.optimize import *
# from numdifftools import Jacobian
from math import sqrt,pi,pow
from numpy import *
from scipy.constants import *
from scipy.signal import *
from numpy import linalg
#from scipy import *
from numpy import fft
# from time import time
import numpy

from EXAFS_analysis.xafs.FT import *
from EXAFS_analysis.xafs.def_sample import *
from EXAFS_analysis.xafs.FEFF_cal import *
from EXAFS_analysis.xafs.ratio_fit import *
from EXAFS_analysis.xafs.artemis_FITplot import *

def interpolation(f, x):
    g = interp1d(f[:, 0], f[:, 1])
    return g(x.tolist())

