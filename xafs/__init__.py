#!/usr/bin/env python

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import glob
from scipy import stats
from scipy.optimize import *
from numdifftools import Jacobian

from xafs.FT import *
from xafs.def_sample import *
from xafs.FEFF_cal import *
from xafs.ratio_fit import *
from xafs.artemis_FITplot import *