#!/usr/bin/env python

from xafs import *



class exafs_func():
    params_dic = {'O': [3, 2.4, 0.002], 'Cd': [1, 4.0, 0.02], 'S': [2.5, 2.52, 0.0031], 'C':[4,3.5,0.01], 'sd': [1], 'e0': [0]}
    bounds_dic = {'S': np.array([(0, 4), (2.4, 2.6), (0.002, 0.007)]),
                  'Cd': np.array([(0, 12), (3.8, 4.3), (0.003, 0.1)]),
                  'O': np.array([(0, 6), (2.1, 2.5), (0.001, 0.02)]),
                  'C': np.array([(0, 8), (2.4, 4), (0.001, 0.2)]),
                  'sd': np.array([(0.5, 1.5)]),
                  'e0': np.array([(0,0)])}

    def __init__(self, sample='bulk', k_min=2, k_max=18, R_CdO=2.335, R_CdS=2.5256, ss_CdS = 0.0031, ss_CdO = 0.0031):


        self.R_CdO = R_CdO
        self.R_CdS = R_CdS
        self.ss_CdS = ss_CdS
        self.ss_CdO = ss_CdO
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
        FEFF_CdO = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/CdO/feff0001.dat',
                                 skip_header=15)

        # CdC path data : Cd -> O -> C bond
        FEFF_CdC = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/CdCO3/feff0002.dat',
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

        self.k_min = k_min
        self.k_max = k_max

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
        plt.subplot(2,1,1)
        for i in filelist:
            data = np.genfromtxt(i)[:, (0, 2)]
            data = self.interpolation(data, self.k)
            plt.plot(self.k, data * self.k ** 2, label='Artemis {:s}'.format(i[59:-2]))
        plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)
        plt.subplot(2,1,2)
        for i in filelist:
            data = np.genfromtxt(i)[:, (0, 2)]
            data = self.interpolation(data, self.k)
            plt.plot(self.k, data * self.k ** 2 - self.experiment * self.k ** 2, label='Artemis {:s}'.format(i[59:-2]))
        plt.legend()
        plt.xlabel('k')
        plt.ylabel('Chi')
        plt.title('Artemis fit and {:s} experiment'.format(self.sample))
        plt.show()
        return


    def EXAFS_model_COS(self, x, N1, R1, ss1, N3, R3, ss3, N4, R4,ss4):
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
        return N1 * x.astype(float) ** 2 \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N2 * x.astype(float) ** 2 \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
               * np.sin(2 * x.astype(float) * R2 + self.phase2.astype(float)) \
               * np.exp(-2 * R2 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss2 * x.astype(float) ** 2)

    def EXAFS_model_COSCd(self, x, N1, R1, ss1, N2, R2, ss2, N3, R3, ss3, N4, R4, ss4):
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
               * np.exp(-2 * ss3 * x.astype(float) ** 2) \
               + N4 * x ** 2 \
               * abs(self.Feff4.astype(float)) / (x.astype(float) * R4 ** 2) \
               * np.sin(2 * x.astype(float) * R4 + self.phase4.astype(float)) \
               * np.exp(-2 * R4 / self.lambda4.astype(float)) \
               * np.exp(-2 * ss4 * x.astype(float) ** 2)

    # def dg(self, param, data, fit_model):
    #     def g_shell(x,y,param1,param2,param3,feff,phase,lambdaX,f):
    #         g1 = np.sum(x ** 2
    #                     * abs(feff.astype(float)) / (x.astype(float) * param2 ** 2)
    #                     * np.sin(2 * x.astype(float) * param2 + phase.astype(float))
    #                     * np.exp(-2 * param2 / lambdaX.astype(float))
    #                     * np.exp(-2 * param3 * x.astype(float) ** 2)
    #                     * (f(x, *param)-y))
    #         g2 = np.sum((-2 * param1 * x ** 2
    #                     * abs(feff.astype(float)) / (x.astype(float) * param2 ** 3)
    #                     * np.sin(2 * x.astype(float) * param2 + phase.astype(float))
    #                     * np.exp(-2 * param2 / lambdaX.astype(float))
    #                     * np.exp(-2 * param3 * x.astype(float) ** 2)
    #                     + 2 * x.astype(float) * param1 * x ** 2
    #                     * abs(feff.astype(float)) / (x.astype(float) * param2 ** 2)
    #                     * np.cos(2 * x.astype(float) * param2 + phase.astype(float))
    #                     * np.exp(-2 * param2 / lambdaX.astype(float))
    #                     * np.exp(-2 * param3 * x.astype(float) ** 2)
    #                     - 2 / self.lambda1.astype(float) * param1 * x ** 2
    #                     * abs(feff.astype(float)) / (x.astype(float) * param2 ** 2)
    #                     * np.sin(2 * x.astype(float) * param2 + phase.astype(float))
    #                     * np.exp(-2 * param2 / lambdaX.astype(float))
    #                     * np.exp(-2 * param3 * x.astype(float) ** 2))
    #                     * (f(x, *param)-y))
    #         g3 = np.sum((-2 * x **2 * param1 * x ** 2
    #                     * abs(feff.astype(float)) / (x.astype(float) * param2 ** 2)
    #                     * np.sin(2 * x.astype(float) * param2 + phase.astype(float))
    #                     * np.exp(-2 * param2 / lambdaX.astype(float))
    #                     * np.exp(-2 * param3 * x.astype(float) ** 2))
    #                     * (f(x, *param)-y))
    #         return g1,g2,g3
    #
    #     x, y = data
    #
    #     if fit_model == 'COS':
    #         N1, R1, ss1, N3, R3, ss3, N4, R4, ss4 = param
    #
    #         g1,g2,g3 = g_shell(x,y,N1,R1,ss1,self.Feff1,self.phase1,self.lambda1,self.EXAFS_model_COS)
    #         g4,g5,g6 = g_shell(x,y,N3,R3,ss3,self.Feff3,self.phase3,self.lambda3,self.EXAFS_model_COS)
    #         g7,g8,g9 = g_shell(x,y,N4,R4,ss4,self.Feff4,self.phase4,self.lambda4,self.EXAFS_model_COS)
    #         return g1,g2,g3,g4,g5,g6,g7,g8,g9
    #     if fit_model == 'COSCd':
    #         N1, R1, ss1, N2, R2, ss2, N3, R3, ss3, N4, R4, ss4 = param
    #         g1,g2,g3 = g_shell(x,y,N1,R1,ss1,self.Feff1,self.phase1,self.lambda1,EXAFS_model_COSCd)
    #         g4,g5,g6 = g_shell(x,y,N2,R2,ss2,self.Feff2,self.phase2,self.lambda2,EXAFS_model_COSCd)
    #         g7,g8,g9 = g_shell(x,y,N3,R3,ss3,self.Feff3,self.phase3,self.lambda3,EXAFS_model_COSCd)
    #         g10,g11,g12 = g_shell(x,y,N4,R4,ss4,self.Feff4,self.phase4,self.lambda4,EXAFS_model_COSCd)
    #         return g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12

    ###### jacobian in non-linear fit
    # def CD(self, f, x0, h):
    #     """Returns the derivative of f at x0 evaluated using the central difference algorithm with step size h."""
    #     return (f(x0 + h/2) - f(x0 - h/2))/np.linalg.norm(h)
    #
    #
    # def jacobian(self,param,data,fit_model):
    #     """Returns the Jacobian matrix of g evaluated at param, given observed data."""
    #     p = np.array(param)
    #     delta = 1e-6
    #     N = len(param)
    #
    #     # Start with an empty matrix of the right size.
    #     jac = np.zeros((N, N))
    #
    #     # We want to calculate df_i/dp_j for all i and j, so need two loops.
    #     for i in range(N):
    #         # Define an appropriate one-dimensional function f_i:
    #         g = self.dg
    #         def g_i(x):
    #             return g(x, data, fit_model)[i]
    #         for j in range(N):
    #             # Set up a step of delta in the appropriate direction:
    #             dv = np.zeros(N)
    #             dv[j] = delta
    #
    #             jac[i,j] = self.CD(g_i, p, dv)
    #
    #     return jac

    def NLcurve_fit(self):

        bounds_OSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['O']))
        params_OSCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['O']
        bounds_OS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O']))
        params_OS = self.params_dic['S'] + self.params_dic['O']
        bounds_SCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd']))
        params_SCd = self.params_dic['S'] + self.params_dic['Cd']
        bounds_COS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['C']))
        params_COS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['C']
        bounds_COSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['C'], self.bounds_dic['Cd']))
        params_COSCd = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['C'] + self.params_dic['Cd']
        plt.subplot(2,1,1)
        # SCd shell
        popt, pcov = curve_fit(self.EXAFS_model_SCd, self.k.astype(float),
                               (self.experiment * self.k ** 2).astype(float),
                               bounds=bounds_SCd.transpose(), p0=params_SCd)
        print('SCd: ',popt)
        fit_y1 = self.EXAFS_model_SCd(self.k, *popt)
        plt.plot(self.k, fit_y1, label='nonlinear curve fit (scipy) with S Cd')
        # OS shell
        popt, pcov = curve_fit(self.EXAFS_model_OS, self.k.astype(float), (self.experiment * self.k ** 2).astype(float),
                               bounds=bounds_OS.transpose(), p0=params_OS)
        print('OS: ',popt)
        fit_y2 = self.EXAFS_model_OS(self.k, *popt)
        plt.plot(self.k, fit_y2, label='nonlinear curve fit (scipy) with O S')
        # OSCd shell
        popt, pcov = curve_fit(self.EXAFS_model_OSCd, self.k.astype(float),
                               (self.experiment * self.k ** 2).astype(float),
                               bounds=bounds_OSCd.transpose(), p0=params_OSCd)
        print('OSCd: ',popt)
        fit_y3 = self.EXAFS_model_OSCd(self.k, *popt)
        plt.plot(self.k, fit_y3, label='nonlinear curve fit (scipy) with O S Cd')
        # COS shell
        popt, pcov = curve_fit(self.EXAFS_model_COS, self.k.astype(float),
                               (self.experiment * self.k ** 2).astype(float),
                               bounds=bounds_COS.transpose(), p0=params_COS)
        print('COS: ',popt)
        fit_y4 = self.EXAFS_model_COS(self.k, *popt)
        plt.plot(self.k, fit_y4, label='nonlinear curve fit (scipy) with C O S')

        # COSCd shell
        popt, pcov = curve_fit(self.EXAFS_model_COSCd, self.k.astype(float),
                               (self.experiment * self.k ** 2).astype(float),
                               bounds=bounds_COSCd.transpose(), p0=params_COSCd)
        print('COSCd: ',popt)
        fit_y5 = self.EXAFS_model_COSCd(self.k, *popt)
        plt.plot(self.k, fit_y3, label='nonlinear curve fit (scipy) with C O S Cd')

        # EXPERIMENT
        plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)

        plt.subplot(2,1,2)
        plt.plot(self.k, fit_y1 - self.experiment * self.k ** 2, label='nonlinear curve fit (scipy) with S Cd')
        plt.plot(self.k, fit_y2 - self.experiment * self.k ** 2, label='nonlinear curve fit (scipy) with O S')
        plt.plot(self.k, fit_y3 - self.experiment * self.k ** 2, label='nonlinear curve fit (scipy) with O S Cd')
        plt.plot(self.k, fit_y4 - self.experiment * self.k ** 2, label='nonlinear curve fit (scipy) with C O S')
        plt.plot(self.k, fit_y5 - self.experiment * self.k ** 2, label='nonlinear curve fit (scipy) with C O S Cd')

        plt.figlegend()
        plt.xlabel('k')
        plt.ylabel('Chi')
        plt.title('nonlinear curve fitting (least square)')
        plt.show()
        return




    def LLE_model_OSCd(self, params):
        model = self.EXAFS_model_OSCd(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** 2, loc=model))
        return LL

    def LLE_model_OS(self, params):
        model = self.EXAFS_model_OS(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** 2, loc=model))
        return LL

    def LLE_model_SCd(self, params):
        model = self.EXAFS_model_SCd(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment.astype(float) * self.k ** 2, loc=model))
        return LL

    def LLE_model_COS(self, params):
        model = self.EXAFS_model_COS(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** 2, loc=model))
        return LL

    def LLE_model_COSCd(self, params):
        model = self.EXAFS_model_COSCd(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** 2, loc=model))
        return LL


    def LLE_fit(self):
        bounds_OSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['O'],
                                 self.bounds_dic['sd']))
        params_OSCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['O'] + self.params_dic['sd']
        bounds_OS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['sd']))
        params_OS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['sd']
        bounds_SCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['sd']))
        params_SCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['sd']
        bounds_COS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['C'],
                                self.bounds_dic['sd']))
        params_COS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['C'] + self.params_dic['sd']
        bounds_COSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['C'],
                                  self.bounds_dic['Cd'], self.bounds_dic['sd']))
        params_COSCd = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['C'] + self.params_dic['Cd'] \
                       + self.params_dic['sd']
        # print(params_SCd)
        plt.figure()
        # SCd shell

        LLE_fitresult1 = minimize(self.LLE_model_SCd, params_SCd, method='trust-constr', bounds=bounds_SCd,
                                  options={'maxiter': 12000})
        print('SCd: ', LLE_fitresult1)
        fit_y1 = self.EXAFS_model_SCd(self.k, *LLE_fitresult1.x[:-1])
        plt.plot(self.k, fit_y1, label='LLE = {:f} for SCd shell'.format(self.LLE_model_SCd(LLE_fitresult1.x)))

        #OS shell
        LLE_fitresult2 = minimize(self.LLE_model_OS, params_OS, method='trust-constr',
                                  bounds=bounds_OS, options={'maxiter': 12000})
        print('OS: ', LLE_fitresult2)
        fit_y2 = self.EXAFS_model_OS(self.k, *LLE_fitresult2.x[:-1])
        plt.plot(self.k, fit_y2, label='LLE = {:f} for OS shell'.format(self.LLE_model_OS(LLE_fitresult2.x)))
        # OSCd shell
        LLE_fitresult3 = minimize(self.LLE_model_OSCd, params_OSCd, method='trust-constr',
                                  bounds=bounds_OSCd, options={'maxiter': 12000})
        print('OSCd: ', LLE_fitresult3)
        fit_y3 = self.EXAFS_model_OSCd(self.k, *LLE_fitresult3.x[:-1])
        plt.plot(self.k, fit_y3, label='LLE = {:f} for OSCd shell'.format(self.LLE_model_OSCd(LLE_fitresult3.x)))

        # jacobian = self.jacobian
        #
        # def jac():
        #     return jacobian(params_COS,[self.k,self.experiment],'COS')

        # COS shell
        LLE_fitresult4 = minimize(self.LLE_model_COS, params_COS, options={'maxiter': 12000},method='trust-constr',
                                   bounds=bounds_COS)
        print('COS: ', LLE_fitresult4)
        fit_y4 = self.EXAFS_model_COS(self.k, *LLE_fitresult4.x[:-1])
        plt.plot(self.k, fit_y4, label='LLE = {:f} for COS shell'.format(self.LLE_model_COS(LLE_fitresult4.x)))

        # COSCd shell
        LLE_fitresult5 = minimize(self.LLE_model_COSCd, params_COSCd, method='trust-constr',
                                  bounds=bounds_COSCd, options={'maxiter': 12000})
        print('COSCd: ', LLE_fitresult5)
        fit_y5 = self.EXAFS_model_COSCd(self.k, *LLE_fitresult5.x[:-1])
        plt.plot(self.k, fit_y5, label='LLE = {:f} for COSCd shell'.format(self.LLE_model_COSCd(LLE_fitresult5.x)))

        # Experiment
        plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=2, label=self.sample)
        plt.title('Min Log Likelihood Estimate (LLE) curve fit')
        plt.legend()
        plt.xlabel('K')
        plt.ylabel('Chi')
        plt.show()

        # difference between fit and experiment
        plt.figure()
        plt.plot(self.k, fit_y1 - self.experiment * self.k ** 2, label='SCd')
        plt.plot(self.k, fit_y2 - self.experiment * self.k ** 2, label='OS')
        plt.plot(self.k, fit_y3 - self.experiment * self.k ** 2, label='OSCd')
        plt.plot(self.k, fit_y4 - self.experiment * self.k ** 2, label='COS')
        plt.plot(self.k, fit_y5 - self.experiment * self.k ** 2, label='COSCd')
        plt.title('difference')
        plt.legend()
        plt.xlabel('K')
        plt.ylabel('Chi')
        plt.show()

        # FT of difference
        plt.figure()
        #
        r,amp = FT_chi(self.k, fit_y1 - self.experiment * self.k ** 2,dx=3)
        plt.plot(r, amp, label='SCd')

        r,amp = FT_chi(self.k, fit_y2 - self.experiment * self.k ** 2,dx=3)
        plt.plot(r, amp, label='OS')

        r,amp = FT_chi(self.k, fit_y3 - self.experiment * self.k ** 2,dx=3)
        plt.plot(r, amp, label='OSCd')

        r,amp = FT_chi(self.k, fit_y4 - self.experiment * self.k ** 2,dx=3)
        plt.plot(r, amp, label='COS')

        r,amp = FT_chi(self.k, fit_y5 - self.experiment * self.k ** 2,dx=3)
        plt.plot(r, amp, label='COSCd')
        plt.xlim([0,5])
        plt.legend()
        plt.title('FT of difference')
        #plt.ylim([0,1])
        plt.xlabel('R')
        plt.ylabel('amp')

        plt.show()
        return fit_y1, fit_y2, fit_y3, fit_y4, fit_y5












