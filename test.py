import numpy as np
from xafs.scattering_path import scattering_path
from xafs.fit import fit

data = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Statistic Analysis/CdS/experiment/CdS_R_Nov17.chik')

CdO = scattering_path('CdO')
CdS = scattering_path('CdS')
CdCd = scattering_path('CdCd')
CdO.cal_exp_params()
CdS.cal_exp_params()
CdCd.cal_exp_params()

R_exp = fit(data[:,0],data[:,1])
R_exp.set_fit_range(2.5,16)
R_exp.set_shells([CdO,CdS,CdCd])
R_exp.set_sample_name('R')
R_exp.fit_loop(loop_num=10)
R_exp.fit_result()
