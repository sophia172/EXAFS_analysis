import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import glob
from scipy import stats
from scipy.optimize import minimize


class exafs_func():
    params_dic = {'O': [5, 2.45, 0.002], 'Cd': [6, 4.1, 0.02], 'S': [3, 2.57, 0.002], 'C':[2.5,3.5,0.01], 'sd': [0.1]}
    bounds_dic = {'S': np.array([(1, 4), (2.3, 2.62), (0.001, 0.005)]), \
                  'Cd': np.array([(1, 12), (3.8, 4.3), (0.01, 0.2)]), \
                  'O': np.array([(0, 6), (2.3, 2.6), (0.001, 0.01)]), \
                  'C': np.array([(0, 6), (2.5, 4), (0.001, 0.02)]), \
                  'sd': np.array([(0.00000001, 0.2)])}

    def __init__(self, sample='bulk', k_min=2, k_max=18):
        self.sample = sample

        if sample == 'bulk':
            exp_col = 2
        elif sample == 'M311':
            exp_col = 3
        else:
            exp_col = 4

        experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt')
        k = experiments[:, 0]
        self.experiment = experiments
        ####
        #
        #     Required data in FEFF analysis
        #

        FEFF_CdS = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/feff0001.dat',
                                 skip_header=15)
        FEFF_CdCd = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/F-43m/feff0002.dat',
                                  skip_header=15)
        # CdO path data
        FEFF_CdO = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/XANES/CdO/feff0001.dat',
                                 skip_header=15)

        # CdC path data : Cd -> O -> C bond
        FEFF_CdC = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/XANES/CdCO3/feff0002.dat',
                                 skip_header=16)

        # feed data files to variable
        Feff1 = FEFF_CdS[:, (0, 2)]
        Feff2 = FEFF_CdCd[:, (0, 2)]
        Feff3 = FEFF_CdO[:, (0, 2)]
        Feff4 = FEFF_CdC[:, (0, 2)]

        phase1 = FEFF_CdS[:, (0, 1)]
        phase1[:, 1] += FEFF_CdS[:, 3]

        phase2 = FEFF_CdCd[:, (0, 1)]
        phase2[:, 1] += FEFF_CdCd[:, 3]

        phase3 = FEFF_CdO[:, (0, 1)]
        phase3[:, 1] += FEFF_CdO[:, 3]

        phase4 = FEFF_CdC[:, (0, 1)]
        phase4[:, 1] += FEFF_CdC[:, 3]

        lambda1 = FEFF_CdS[:, (0, 5)]
        lambda2 = FEFF_CdCd[:, (0, 5)]
        lambda3 = FEFF_CdO[:, (0, 5)]
        lambda4 = FEFF_CdC[:, (0, 5)]

        S01 = FEFF_CdS[:, (0, 4)]
        S02 = FEFF_CdCd[:, (0, 4)]
        S03 = FEFF_CdO[:, (0, 4)]
        S04 = FEFF_CdC[:, (0, 4)]

        k_min = (np.abs(self.experiment[:, 0] - k_min)).argmin()
        k_max = (np.abs(self.experiment[:, 0] - k_max)).argmin()

        self.k = self.experiment[k_min:k_max, 0]
        self.experiment = self.experiment[k_min:k_max, exp_col]

        # print('k = ',k)
        # interpolate given parameter files
        self.Feff1 = self.interpolation(Feff1, self.k)
        self.Feff2 = self.interpolation(Feff2, self.k)
        self.Feff3 = self.interpolation(Feff3, self.k)
        self.Feff4 = self.interpolation(Feff4, self.k)
        self.lambda1 = self.interpolation(lambda1, self.k)
        self.lambda2 = self.interpolation(lambda2, self.k)
        self.lambda3 = self.interpolation(lambda3, self.k)
        self.lambda4 = self.interpolation(lambda4, self.k)
        self.phase1 = self.interpolation(phase1, self.k)
        self.phase2 = self.interpolation(phase2, self.k)
        self.phase3 = self.interpolation(phase3, self.k)
        self.phase4 = self.interpolation(phase4, self.k)
        self.S01 = self.interpolation(S01, self.k)
        self.S02 = self.interpolation(S02, self.k)
        self.S03 = self.interpolation(S03, self.k)
        self.S04 = self.interpolation(S04, self.k)
        return

    def interpolation(self, f, x):
        g = interp1d(f[:, 0], f[:, 1])
        return g(x.tolist())

    ####    Plot artemis fit

    def plot_artemis(self):
        filelist = glob.glob('/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/{:s}*.k'.format(self.sample))
        for i in filelist:
            data = np.genfromtxt(i)[:, (0, 2)]
            data = self.interpolation(data, self.k)
            plt.plot(self.k, data * self.k ** 2, label='Artemis {:s}'.format(i[59:-2]))
        plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)
        plt.legend()
        plt.xlabel('k')
        plt.ylabel('Chi')
        plt.title('Artemis fit and experiment')
        plt.show()
        return

    def plot_artemis_diff(self):
        filelist = glob.glob('/Users/Sophia/ownCloud/PhD/Statistic Analysis/demeter/{:s}*.k'.format(self.sample))
        for i in filelist:
            data = np.genfromtxt(i)[:, (0, 2)]
            data = self.interpolation(data, self.k)
            plt.plot(self.k, data * self.k ** 2 - self.experiment * self.k ** 2, label='Artemis {:s}'.format(i[59:-2]))
        plt.legend()
        plt.xlabel('k')
        plt.ylabel('Chi difference')
        plt.title('difference between Artemis fit and {:s} experiment'.format(self.sample))
        plt.show()
        return

    def plot_xmus(self, folder):
        filelist = glob.glob('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/Cd/{:s}/*xmu.dat'.format(folder))
        for file in filelist:
            data = np.genfromtxt(file)
            plt.plot(data[:, 2], data[:, 5] * data[:, 2] ** 2)
        plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)
        plt.legend()
        plt.xlabel('k')
        plt.ylabel('Chi')
        plt.title('FEFF InP model (SCd shell)')
        plt.show()
        return
    def EXAFS_model_COS(self,x, N1, R1, ss1, N3, R3, ss3, N4, R4,ss4):
        return N1 * x ** 2 \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N3 * x ** 2 \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
               * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
               * np.exp(-2 * R3 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss3 * x.astype(float) ** 2) \
               + N4 * x ** 2 \
               * abs(self.Feff4.astype(float)) / (x.astype(float) * R4 ** 2) \
               * np.sin(2 * x.astype(float) * R4 + self.phase4.astype(float)) \
               * np.exp(-2 * R4 / self.lambda4.astype(float)) \
               * np.exp(-2 * ss4 * x.astype(float) ** 2)


    def EXAFS_model_OS(self, x, N1, R1, ss1, N3, R3, ss3):
        return N1 * x ** 2 \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N3 * x ** 2 \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
               * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
               * np.exp(-2 * R3 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss3 * x.astype(float) ** 2)

    def EXAFS_model_OSCd(self, x, N1, R1, ss1, N2, R2, ss2, N3, R3, ss3):
        return N1 * x ** 2 \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N2 * x ** 2 \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
               * np.sin(2 * x.astype(float) * R2 + self.phase2.astype(float)) \
               * np.exp(-2 * R2 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss2 * x.astype(float) ** 2) \
               + N3 * x ** 2 \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
               * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
               * np.exp(-2 * R3 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss3 * x.astype(float) ** 2)

    def EXAFS_model_SCd(self, x, N1, R1, ss1, N2, R2, ss2):
        return N1 * x ** 2 \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N2 * x ** 2 \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
               * np.sin(2 * x.astype(float) * R2 + self.phase2.astype(float)) \
               * np.exp(-2 * R2 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss2 * x.astype(float) ** 2)

    def NLcurve_fit(self):

        bounds_OSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['O']))
        params_OSCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['O']
        bounds_OS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O']))
        params_OS = self.params_dic['S'] + self.params_dic['O']
        bounds_SCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd']))
        params_SCd = self.params_dic['S'] + self.params_dic['Cd']
        bounds_COS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['C']))
        params_COS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['C']
        # SCd shell
        popt, pcov = curve_fit(self.EXAFS_model_SCd, self.k.astype(float),
                               (self.experiment * self.k ** 2).astype(float), \
                               bounds=bounds_SCd.transpose(), p0=params_SCd)
        print('SCd: ',pcov)
        plt.plot(self.k, self.EXAFS_model_SCd(self.k, *popt), label='nonlinear curve fit (scipy) with S Cd')
        # OS shell
        popt, pcov = curve_fit(self.EXAFS_model_OS, self.k.astype(float), (self.experiment * self.k ** 2).astype(float), \
                               bounds=bounds_OS.transpose(), p0=params_OS)
        print('OS: ',pcov)
        plt.plot(self.k, self.EXAFS_model_OS(self.k, *popt), label='nonlinear curve fit (scipy) with O S')
        # OSCd shell
        popt, pcov = curve_fit(self.EXAFS_model_OSCd, self.k.astype(float),
                               (self.experiment * self.k ** 2).astype(float), \
                               bounds=bounds_OSCd.transpose(), p0=params_OSCd)
        print('OSCd: ',pcov)
        plt.plot(self.k, self.EXAFS_model_OSCd(self.k, *popt), label='nonlinear curve fit (scipy) with O S Cd')
        # COS shell
        popt, pcov = curve_fit(self.EXAFS_model_COS, self.k.astype(float),
                               (self.experiment * self.k ** 2).astype(float), \
                               bounds=bounds_COS.transpose(), p0=params_COS)
        print('COS: ',pcov)
        plt.plot(self.k, self.EXAFS_model_COS(self.k, *popt), label='nonlinear curve fit (scipy) with C O S')
        # EXPERIMENT
        plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)

        plt.legend()
        plt.xlabel('k')
        plt.ylabel('Chi')
        plt.title('nonlinear curve fitting (least square)')
        plt.show()
        return

    def LLE_model_OSCd(self, params):
        sd = params[-1]
        model = self.EXAFS_model_OSCd(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** 2, loc=model, scale=sd))
        return LL

    def LLE_model_OS(self, params):
        sd = params[-1]
        model = self.EXAFS_model_OS(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** 2, loc=model, scale=sd))
        return LL

    def LLE_model_SCd(self, params):
        sd = params[-1]
        model = self.EXAFS_model_SCd(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** 2, loc=model, scale=sd))
        return LL

    def LLE_model_COS(self, params):
        sd = params[-1]
        model = self.EXAFS_model_COS(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** 2, loc=model, scale=sd))
        return LL

    def LLE_fit(self):
        bounds_OSCd = np.vstack(
            (self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['O'], self.bounds_dic['sd']))
        params_OSCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['O'] + self.params_dic['sd']
        bounds_OS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['sd']))
        params_OS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['sd']
        bounds_SCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['sd']))
        params_SCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['sd']
        bounds_COS = np.vstack(
            (self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['C'], self.bounds_dic['sd']))
        params_COS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['C'] + self.params_dic['sd']
        # SCd shell
        LLE_fitresult1 = minimize(self.LLE_model_SCd, params_SCd, bounds=bounds_SCd, \
                                  method='TNC', options={'maxiter': 8000})
        print('SCd: ', LLE_fitresult1)
        plt.plot(self.k, self.EXAFS_model_SCd(self.k, *LLE_fitresult1.x[:-1]), \
                 label='LLE = {:f} for SCd shell'.format(self.LLE_model_SCd(LLE_fitresult1.x)))
        # OS shell
        LLE_fitresult2 = minimize(self.LLE_model_OS, params_OS, bounds=bounds_OS, \
                                  method='TNC', options={'maxiter': 8000})
        print('OS: ', LLE_fitresult2)
        plt.plot(self.k, self.EXAFS_model_OS(self.k, *LLE_fitresult2.x[:-1]), \
                 label='LLE = {:f} for OS shell'.format(self.LLE_model_OS(LLE_fitresult2.x)))
        # OSCd shell
        LLE_fitresult3 = minimize(self.LLE_model_OSCd, params_OSCd, bounds=bounds_OSCd, \
                                  method='TNC', options={'maxiter': 8000})
        print('OSCd: ', LLE_fitresult3)
        plt.plot(self.k, self.EXAFS_model_OSCd(self.k, *LLE_fitresult3.x[:-1]), \
                 label='LLE = {:f} for OSCd shell'.format(self.LLE_model_OSCd(LLE_fitresult3.x)))
        # COS shell
        LLE_fitresult4 = minimize(self.LLE_model_COS, params_COS, bounds=bounds_COS, \
                                  method='TNC', options={'maxiter': 8000})
        print('COS: ', LLE_fitresult4)
        plt.plot(self.k, self.EXAFS_model_COS(self.k, *LLE_fitresult4.x[:-1]), \
                 label='LLE = {:f} for COS shell'.format(self.LLE_model_COS(LLE_fitresult4.x)))
        # Experiment
        plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)

        plt.legend()
        plt.xlabel('k')
        plt.ylabel('Chi')
        plt.title('Min Log Likelihood Estimate (LLE) curve fit')
        plt.show()
        return
