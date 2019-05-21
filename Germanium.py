import numpy as np
from scipy.interpolate import interp1d
from xafs import *
from scipy.optimize import curve_fit

def interpolation(f, x):
    g = interp1d(f[:, 0], f[:, 1])
    return g(x.tolist())

experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi.dat',skip_header=5)
k = experiments[:,0]
exp = experiments[:,1]
weight=2

k_min = (np.abs(k - 3)).argmin()
k_max = (np.abs(k - 19)).argmin()
k = k[k_min:k_max]
exp = exp[k_min:k_max]

FEFF_1 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/feff0001.dat',
                         skip_header=15)
FEFF_2 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/feff0002.dat',
                       skip_header=15)
FEFF_3 = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/Ge/feff0005.dat',
                       skip_header=15)

Feff1 = FEFF_1[:, (0, 2)]
Feff2 = FEFF_2[:, (0, 2)]
Feff3 = FEFF_3[:, (0, 2)]


phase1 = FEFF_1[:, (0, 1)]
phase1[:, 1] += FEFF_1[:, 3]

phase2 = FEFF_2[:, (0, 1)]
phase2[:, 1] += FEFF_2[:, 3]

phase3 = FEFF_3[:, (0, 1)]
phase3[:, 1] += FEFF_3[:, 3]


lambda1 = FEFF_1[:, (0, 5)]
lambda2 = FEFF_2[:, (0, 5)]
lambda3 = FEFF_3[:, (0, 5)]

Feff1 =  interpolation(Feff1,  k)
Feff2 =  interpolation(Feff2,  k)
Feff3 =  interpolation(Feff3,  k)
lambda1 =  interpolation(lambda1,  k)
lambda2 =  interpolation(lambda2,  k)
lambda3 =  interpolation(lambda3,  k)
phase1 =  interpolation(phase1,  k)
phase2 =  interpolation(phase2,  k)
phase3 =  interpolation(phase3,  k)
#
# plt.plot(k,phase3)
# plt.show()
print('R1              ss1               R2             ss2               R3                  ss3')
def EXAFS_model(x,N1,R1,ss1,N2,R2,ss2,N3,R3,ss3,e):
    # print(R1,ss1,R2,ss2,R3,ss3)
    return e + N1 * x **weight\
           * abs( Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
           * np.sin(2 * x.astype(float) * R1 +  phase1.astype(float)) \
           * np.exp(-2 * R1 /  lambda1.astype(float)) \
           * np.exp(-2 * ss1 * x.astype(float) ** 2) \
           + N2 * x **weight\
           * abs( Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
           * np.sin(2 * x.astype(float) * R2 +  phase2.astype(float)) \
           * np.exp(-2 * R2 /  lambda2.astype(float)) \
           * np.exp(-2 * ss2 * x.astype(float) ** 2)\
           + N3 * x **weight\
           * abs( Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
           * np.sin(2 * x.astype(float) * R3 +  phase3.astype(float)) \
           * np.exp(-2 * R3 /  lambda3.astype(float)) \
           * np.exp(-2 * ss3 * x.astype(float) ** 2)


params_dic = {'Ge1': [6, 2.49, 0.0058], 'Ge2': [12, 3.77, 0.009], 'Ge3': [12, 4.055, 0.0063], 'e':[0], 'sd': [1]}
bounds_dic = {'Ge1': np.array([(6, 6), (2.4, 2.6), (0.0001, 0.02)]),
              'Ge2': np.array([(12,12), (3.2, 4.3), (0.0001, 0.02)]),
              'Ge3': np.array([(12,12), (3.5, 4.9), (0.0001, 0.02)]),
              'e': np.array([(-10,10)]),
              'sd': np.array([(0,1)])}


params = params_dic['Ge1'] + params_dic['Ge2'] + params_dic['Ge3'] + params_dic['e'] + params_dic['sd']
bounds = np.vstack((bounds_dic['Ge1'], bounds_dic['Ge2'], bounds_dic['Ge3'], bounds_dic['e'], bounds_dic['sd']))

popt, pcov = curve_fit(EXAFS_model, k.astype(float),(exp * k **2).astype(float), p0=params[:-1])
fit = EXAFS_model(k, *popt)
print(*popt)

def LLE_model(params):
    model = EXAFS_model(k, *params[:-1])
    LL = -np.sum(stats.norm.logpdf(exp * k**weight, loc=model))
    return LL
LLE_fitresult = minimize(LLE_model, params, method='TNC', bounds=bounds,options={'maxiter': 12000})
fit_LLE = EXAFS_model(k, *LLE_fitresult.x[:-1])
print(*LLE_fitresult.x[:-1])

plt.figure()
plt.plot(k,exp *k**weight,'k', linewidth = 1, label = 'experiment')
plt.plot(k, fit, linewidth=1, label='nonlinear curve fit (scipy) of Germanium')
plt.plot(k,fit_LLE,linewidth=1, label='LLE fit (scipy) of Germanium')
plt.legend()
plt.show()

plt.figure()
r,amp,real,imag = calcfft(k,exp*k**3,kmin=2,kmax=20)
plt.plot(r,amp,'k', linewidth = 1, label = 'experiment')
r,amp,real,imag = calcfft(k,fit*k,kmin=2,kmax=20)
plt.plot(r,amp, linewidth=1, label='nonlinear curve fit (scipy) of Germanium')
r,amp,real,imag = calcfft(k,fit_LLE*k,kmin=2,kmax=20)
plt.plot(r,amp,linewidth=1, label='LLE fit (scipy) of Germanium')
plt.legend()
plt.xlim([0,5])
plt.show()