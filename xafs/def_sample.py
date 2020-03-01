#!/usr/bin/env python

from EXAFS_analysis.xafs import *



class exafs_func():
    linewidth = 0.5
    weight = 2
    params_dic = {'O': [7.8, 2.4, 0.02528],'O2': [7.8, 2.4, 0.02528], 'Cd': [3, 3.8, 0.06], 'Cd2': [3, 3.4, 0.07],'S': [4, 2.5165, 0.006], 'C':[2,1.5,0.01], 'sd': [1], 'e0': [0]}
    bounds_dic = {'S': np.array([(1, 4), (2.4, 2.6), (0.003, 0.01)]),
                  'Cd': np.array([(1, 12), (3.5, 4.5), (0.01, 0.1)]),
                  'Cd2': np.array([(3, 12), (3.8, 4.5), (0.005, 0.1)]),
                  'O': np.array([(1, 12), (2.1, 2.5), (0.001, 0.05)]),
                  'O2': np.array([(1, 12), (3, 5), (0.001, 0.5)]),
                  'C': np.array([(0, 6), (1, 4), (0.001, 0.2)]),
                  'sd': np.array([(0.5, 1.5)]),
                  'e0': np.array([(-20,20)])}


    def __init__(self, sample='bulk', k_min=2, k_max=18, R_CdO=2.280, R_CdS=2.516, ss_CdS = 0.0047, ss_CdO = 0.0090,fitted_curve='athena'):

        self.fitted_curve = fitted_curve
        self.R_CdO = R_CdO
        self.R_CdS = R_CdS
        self.ss_CdS = ss_CdS
        self.ss_CdO = ss_CdO
        self.sample = sample



        if fitted_curve == 'athena':
            if sample == 'bulk':
                experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_bulk_Aug18_CdK_90K_tran.chiq')
                self.k_original =experiments[1:,0]
                self.experiment_original = experiments[1:,1]/self.k_original**2
            elif sample == 'M311':
                experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt')
                self.k_original = experiments[:, 0]
                self.experiment_original = experiments[:,2]
            else:
                experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_XAFS_CdK_chi_2018.txt')
                self.k_original = experiments[:, 0]
                self.experiment_original = experiments[:,4]
            k_min = (np.abs(self.k_original - k_min)).argmin()
            k_max = (np.abs(self.k_original - k_max)).argmin()
            self.k_max = k_max
            self.k = self.k_original[k_min:k_max]
            self.experiment = self.experiment_original[k_min:k_max]
            # print(self.k,self.experiment)
            # plt.plot(self.k,self.experiment)
            # plt.show()
        ######pyspline chi result######
        #
        ##
        elif fitted_curve == 'pyspline_Andrei':
            # my pyspline
            # experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_{:s}_Aug18_pyspline.dat'.format(sample),skip_header=1)
            # Andrei's pyspline
            experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/andrei_pyspline/CdS_MSC_{:s}_90K_py.dat'.format(sample),skip_header=6)

            self.k_original = experiments[:,2]
            self.experiment_original = experiments[:,7]
            k_min = (np.abs(self.k_original - k_min)).argmin()
            k_max = (np.abs(self.k_original - k_max)).argmin()
            self.k_max = k_max
            self.k = self.k_original[k_min:k_max]
            self.experiment = self.experiment_original[k_min:k_max]

            # print(self.k,self.experiment)
            # plt.plot(self.k,self.experiment)
            # plt.show()

        elif fitted_curve == 'pyspline_Ying':
            # my pyspline
            experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/CdS_{:s}_Aug18_pyspline.dat'.format(sample),skip_header=1)
            # Andrei's pyspline
            # experiments = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Experiment_Analysis/CdS/andrei_pyspline/CdS_MSC_{:s}_90K_py.dat'.format(sample),skip_header=6)

            self.k_original = experiments[:,2]
            self.experiment_original = experiments[:,7]
            k_min = (np.abs(self.k_original - k_min)).argmin()
            k_max = (np.abs(self.k_original - k_max)).argmin()
            self.k_max = k_max
            self.k = self.k_original[k_min:k_max]
            self.experiment = self.experiment_original[k_min:k_max]

            # print(self.k,self.experiment)
            # plt.plot(self.k,self.experiment)
            # plt.show()
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
        FEFF_CdC = np.genfromtxt('/Users/Sophia/ownCloud/PhD/Simulation/CdS/FEFF/EXAFS/CdC/feff0001.dat',
                                 skip_header=15)

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




        # print('k = ',k)
        # interpolate given parameter files, output has only one column

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
        plt.plot(self.k, self.experiment * self.k ** 2, 'k.', markersize=0.5, label=self.sample)
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

    #
    # def EXAFS_model_COS(self, x, N1, R1, ss1, N3, R3, ss3, N4, R4,ss4):
    #     return N1 * x ** self.weight \
    #            * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
    #            * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
    #            * np.exp(-2 * R1 / self.lambda1.astype(float)) \
    #            * np.exp(-2 * ss1 * x.astype(float) ** 2) \
    #            + N3 * x ** self.weight \
    #            * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
    #            * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
    #            * np.exp(-2 * R3 / self.lambda3.astype(float)) \
    #            * np.exp(-2 * ss3 * x.astype(float) ** 2) \
    #            + N4 * x ** self.weight \
    #            * abs(self.Feff4.astype(float)) / (x.astype(float) * R4 ** 2) \
    #            * np.sin(2 * x.astype(float) * R4 + self.phase4.astype(float)) \
    #            * np.exp(-2 * R4 / self.lambda4.astype(float)) \
    #            * np.exp(-2 * ss4 * x.astype(float) ** 2)


    def EXAFS_model_O(self,x,N,R,ss):
        return N * x ** self.weight \
                * abs(self.Feff3.astype(float)) / (x.astype(float) * R ** 2) \
                * np.sin(2 * x.astype(float) * R + self.phase3.astype(float)) \
                * np.exp(-2 * R / self.lambda3.astype(float)) \
                * np.exp(-2 * ss * x.astype(float) ** 2)
    def EXAFS_model_S(self,x,N,R,ss):
            return N * x ** self.weight \
                   * abs(self.Feff1.astype(float)) / (x.astype(float) * R ** 2) \
                   * np.sin(2 * x.astype(float) * R + self.phase1.astype(float)) \
                   * np.exp(-2 * R / self.lambda1.astype(float)) \
                   * np.exp(-2 * ss * x.astype(float) ** 2)
    def EXAFS_model_Cd(self,x,N,R,ss):
        return N * x ** self.weight \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R ** 2) \
               * np.sin(2 * x.astype(float) * R + self.phase2.astype(float)) \
               * np.exp(-2 * R / self.lambda2.astype(float)) \
               * np.exp(-2 * ss * x.astype(float) ** 2)




    def EXAFS_model_OS(self, x, N1, R1, ss1, N3, R3, ss3):
        return N1 * x ** self.weight \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N3 * x ** self.weight \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
               * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
               * np.exp(-2 * R3 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss3 * x.astype(float) ** 2)

    def EXAFS_model_OSCd(self, x, N1, R1, ss1, N2, R2, ss2, N3, R3, ss3):
        return N1 * x ** self.weight \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N2 * x ** self.weight \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
               * np.sin(2 * x.astype(float) * R2 + self.phase2.astype(float)) \
               * np.exp(-2 * R2 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss2 * x.astype(float) ** 2) \
               + N3 * x ** self.weight \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
               * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
               * np.exp(-2 * R3 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss3 * x.astype(float) ** 2)

    def EXAFS_model_OSCd2(self, x, N1, R1, ss1, N2, R2, ss2, N22, R22, ss22, N3, R3, ss3):
        return N1 * x ** self.weight \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N2 * x ** self.weight \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
               * np.sin(2 * x.astype(float) * R2 + self.phase2.astype(float)) \
               * np.exp(-2 * R2 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss2 * x.astype(float) ** 2) \
               + N22 * x ** self.weight \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R22 ** 2) \
               * np.sin(2 * x.astype(float) * R22 + self.phase2.astype(float)) \
               * np.exp(-2 * R22 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss22 * x.astype(float) ** 2) \
               + N3 * x ** self.weight \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
               * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
               * np.exp(-2 * R3 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss3 * x.astype(float) ** 2)


    def EXAFS_model_O2S(self, x, N1, R1, ss1, N3, R3, ss3, N33, R33, ss33):
        return N1 * x ** self.weight \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N3 * x ** self.weight \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
               * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
               * np.exp(-2 * R3 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss33 * x.astype(float) ** 2) \
               + N33 * x ** self.weight \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R33 ** 2) \
               * np.sin(2 * x.astype(float) * R33 + self.phase3.astype(float)) \
               * np.exp(-2 * R33 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss3 * x.astype(float) ** 2)

    def EXAFS_model_O2SCd2(self, x, N1, R1, ss1, N2, R2, ss2, N22, R22, ss22, N3, R3, ss3, N33, R33, ss33):
        return N1 * x ** self.weight \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N2 * x ** self.weight \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
               * np.sin(2 * x.astype(float) * R2 + self.phase2.astype(float)) \
               * np.exp(-2 * R2 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss2 * x.astype(float) ** 2) \
               + N22 * x ** self.weight \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R22 ** 2) \
               * np.sin(2 * x.astype(float) * R22 + self.phase2.astype(float)) \
               * np.exp(-2 * R22 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss22 * x.astype(float) ** 2) \
               + N3 * x ** self.weight \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
               * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
               * np.exp(-2 * R3 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss3 * x.astype(float) ** 2) \
               + N33 * x ** self.weight \
               * abs(self.Feff3.astype(float)) / (x.astype(float) * R33 ** 2) \
               * np.sin(2 * x.astype(float) * R33 + self.phase3.astype(float)) \
               * np.exp(-2 * R33 / self.lambda3.astype(float)) \
               * np.exp(-2 * ss33 * x.astype(float) ** 2)



    def EXAFS_model_SCd(self, x, N1, R1, ss1, N2, R2, ss2):
        return N1 * x.astype(float) ** self.weight \
               * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
               * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
               * np.exp(-2 * R1 / self.lambda1.astype(float)) \
               * np.exp(-2 * ss1 * x.astype(float) ** 2) \
               + N2 * x.astype(float) ** self.weight \
               * abs(self.Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
               * np.sin(2 * x.astype(float) * R2 + self.phase2.astype(float)) \
               * np.exp(-2 * R2 / self.lambda2.astype(float)) \
               * np.exp(-2 * ss2 * x.astype(float) ** 2)
    #
    # def EXAFS_model_COSCd(self, x, N1, R1, ss1, N2, R2, ss2, N3, R3, ss3, N4, R4, ss4):
    #     return N1 * x ** self.weight \
    #            * abs(self.Feff1.astype(float)) / (x.astype(float) * R1 ** 2) \
    #            * np.sin(2 * x.astype(float) * R1 + self.phase1.astype(float)) \
    #            * np.exp(-2 * R1 / self.lambda1.astype(float)) \
    #            * np.exp(-2 * ss1 * x.astype(float) ** 2) \
    #            + N2 * x ** self.weight \
    #            * abs(self.Feff2.astype(float)) / (x.astype(float) * R2 ** 2) \
    #            * np.sin(2 * x.astype(float) * R2 + self.phase2.astype(float)) \
    #            * np.exp(-2 * R2 / self.lambda2.astype(float)) \
    #            * np.exp(-2 * ss2 * x.astype(float) ** 2) \
    #            + N3 * x ** self.weight \
    #            * abs(self.Feff3.astype(float)) / (x.astype(float) * R3 ** 2) \
    #            * np.sin(2 * x.astype(float) * R3 + self.phase3.astype(float)) \
    #            * np.exp(-2 * R3 / self.lambda3.astype(float)) \
    #            * np.exp(-2 * ss3 * x.astype(float) ** 2) \
    #            + N4 * x ** self.weight \
    #            * abs(self.Feff4.astype(float)) / (x.astype(float) * R4 ** 2) \
    #            * np.sin(2 * x.astype(float) * R4 + self.phase4.astype(float)) \
    #            * np.exp(-2 * R4 / self.lambda4.astype(float)) \
    #            * np.exp(-2 * ss4 * x.astype(float) ** 2)

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

        bounds_OSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['O'], self.bounds_dic['e0']))
        params_OSCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['O'] + self.params_dic['e0']
        bounds_OSCd2 = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['Cd'], self.bounds_dic['O'], self.bounds_dic['e0']))
        params_OSCd2 = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['Cd'] + self.params_dic['O'] + self.params_dic['e0']
        bounds_O2S = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['O2'], self.bounds_dic['e0']))
        params_O2S = self.params_dic['S']  + self.params_dic['O'] + self.params_dic['O2'] + self.params_dic['e0']
        bounds_O2SCd2 = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['Cd'], self.bounds_dic['O'], self.bounds_dic['O'], self.bounds_dic['e0']))
        params_O2SCd2 = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['Cd'] + self.params_dic['O'] + self.params_dic['O'] + self.params_dic['e0']

        bounds_OS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['e0']))
        params_OS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['e0']
        bounds_SCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['e0']))
        params_SCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['e0']
        bounds_COS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['C'], self.bounds_dic['e0']))
        params_COS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['C'] + self.params_dic['e0']
        bounds_COSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['O'], self.bounds_dic['C'], self.bounds_dic['e0']))
        params_COSCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['O'] + self.params_dic['C'] + self.params_dic['e0']
        plt.subplot(2,1,1)
        # SCd shell
        popt, pcov = curve_fit(self.EXAFS_model_SCd, self.k.astype(float),
                               (self.experiment * self.k ** self.weight).astype(float),
                               bounds=bounds_SCd.transpose(), p0=params_SCd)
        print('SCd: ',popt)
        fit_y1 = self.EXAFS_model_SCd(self.k, *popt)
        LLE1_NL = self.LLE_model_SCd((*popt,1))
        plt.plot(self.k, fit_y1, linewidth=self.linewidth, label='nonlinear curve fit (scipy) with S Cd')

        # OS shell
        popt, pcov = curve_fit(self.EXAFS_model_OS, self.k.astype(float), (self.experiment * self.k ** self.weight).astype(float),
                               bounds=bounds_OS.transpose(), p0=params_OS)
        print('OS: ',popt)
        fit_y2 = self.EXAFS_model_OS(self.k, *popt)
        LLE2_NL = self.LLE_model_OS((*popt,1))
        plt.plot(self.k, fit_y2, linewidth=self.linewidth, label='nonlinear curve fit (scipy) with O S')

        # OSCd shell
        popt, pcov = curve_fit(self.EXAFS_model_OSCd, self.k.astype(float),
                               (self.experiment * self.k ** self.weight).astype(float),
                               bounds=bounds_OSCd.transpose(), p0=params_OSCd)
        print('OSCd: ',popt)
        fit_y3 = self.EXAFS_model_OSCd(self.k, *popt)
        LLE3_NL = self.LLE_model_OSCd((*popt,1))
        plt.plot(self.k, fit_y3, linewidth=self.linewidth, label='nonlinear curve fit (scipy) with O S Cd')




        # # COS shell
        # popt, pcov = curve_fit(self.EXAFS_model_COS, self.k.astype(float),
        #                        (self.experiment * self.k ** self.weight).astype(float),
        #                        bounds=bounds_COS.transpose(), p0=params_COS)
        # print('COS: ',popt)
        # fit_y4 = self.EXAFS_model_COS(self.k, *popt)
        # LLE4_NL = self.LLE_model_COS((*popt,1))
        # plt.plot(self.k, fit_y4, linewidth=self.linewidth, label='nonlinear curve fit (scipy) with C O S')
        #
        # # COSCd shell
        # popt, pcov = curve_fit(self.EXAFS_model_COSCd, self.k.astype(float),
        #                        (self.experiment * self.k ** self.weight).astype(float),
        #                        bounds=bounds_COSCd.transpose(), p0=params_COSCd)
        # print('COSCd: ',popt)
        # fit_y5 = self.EXAFS_model_COSCd(self.k, *popt)
        # LLE5_NL = self.LLE_model_COSCd((*popt,1))
        # plt.plot(self.k, fit_y5, linewidth=self.linewidth, label='nonlinear curve fit (scipy) with C O S Cd')


        # OSCdCd shell
        popt, pcov = curve_fit(self.EXAFS_model_OSCd2, self.k.astype(float),
                               (self.experiment * self.k ** self.weight).astype(float),
                               bounds=bounds_OSCd2.transpose(),  p0=params_OSCd2)
        print('OSCdCd: ',popt)
        fit_y6 = self.EXAFS_model_OSCd2(self.k, *popt)
        LLE6_NL = self.LLE_model_OSCd2((*popt,1))
        plt.plot(self.k, fit_y6, linewidth=self.linewidth, label='nonlinear curve fit (scipy) with O S Cd Cd')

        # OOSCd shell
        popt, pcov = curve_fit(self.EXAFS_model_O2SCd, self.k.astype(float),
                               (self.experiment * self.k ** self.weight).astype(float),
                               bounds=bounds_O2SCd.transpose(), p0=params_O2SCd)
        print('OOSCd: ',popt)
        fit_y7 = self.EXAFS_model_O2SCd(self.k, *popt)
        LLE7_NL = self.LLE_model_O2SCd((*popt,1))
        plt.plot(self.k, fit_y7, linewidth=self.linewidth, label='nonlinear curve fit (scipy) with O O S Cd')

        # OOSCdCd shell
        popt, pcov = curve_fit(self.EXAFS_model_O2SCd2, self.k.astype(float),
                               (self.experiment * self.k ** self.weight).astype(float),
                               bounds=bounds_O2SCd2.transpose(), p0=params_O2SCd2)
        print('OOSCd: ',popt)
        fit_y8 = self.EXAFS_model_O2SCd2(self.k, *popt)
        LLE8_NL = self.LLE_model_O2SCd2((*popt,1))
        plt.plot(self.k, fit_y8, linewidth=self.linewidth, label='nonlinear curve fit (scipy) with O O S Cd')


        # EXPERIMENT
        plt.plot(self.k, self.experiment * self.k ** self.weight, 'k.', markersize=0.5, label=self.sample)
        plt.xlabel('k')
        plt.ylabel('Chi')
        plt.legend(fontsize=5)
        plt.title('nonlinear curve fitting (least square)')

        plt.subplot(2,1,2)
        plt.plot(self.k, fit_y1 - self.experiment * self.k ** self.weight, label='LLE = {:f} for SCd'.format(LLE1_NL))
        plt.plot(self.k, fit_y2 - self.experiment * self.k ** self.weight, label='LLE = {:f} for OS'.format(LLE2_NL))
        plt.plot(self.k, fit_y3 - self.experiment * self.k ** self.weight, label='LLE = {:f} for OSCd'.format(LLE3_NL))
        # plt.plot(self.k, fit_y4 - self.experiment * self.k ** self.weight, label='LLE = {:f} for COS'.format(LLE4_NL))
        # plt.plot(self.k, fit_y5 - self.experiment * self.k ** self.weight, label='LLE = {:f} for COSCd'.format(LLE5_NL))
        plt.plot(self.k, fit_y6 - self.experiment * self.k ** self.weight, label='LLE = {:f} for OSCdCd'.format(LLE6_NL))
        plt.plot(self.k, fit_y7 - self.experiment * self.k ** self.weight, label='LLE = {:f} for OOSCd'.format(LLE7_NL))
        plt.plot(self.k, fit_y8 - self.experiment * self.k ** self.weight, label='LLE = {:f} for OOSCdCd'.format(LLE8_NL))

        plt.xlabel('k')
        plt.ylabel('Chi')
        plt.title('fit difference')
        plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/fitNL.pdf',format='pdf')
        plt.show()

        # Fourier Transform of difference
        plt.figure()
        r,amp,real,imag = calcfft(self.k, fit_y1 - self.experiment * self.k **3 )
        plt.plot(r, amp,linewidth=self.linewidth, label='SCd')

        r,amp,real,imag = calcfft(self.k, fit_y2 - self.experiment * self.k **3 )
        plt.plot(r, amp, linewidth=self.linewidth, label='OS')

        r,amp,real,imag = calcfft(self.k, fit_y3 - self.experiment * self.k **3 )
        plt.plot(r, amp, linewidth=self.linewidth, label='OSCd')

        # r,amp = calcfft(self.k, fit_y4 - self.experiment * self.k ** self.weight )
        # plt.plot(r, amp, linewidth=self.linewidth, label='COS')
        #
        # r,amp = calcfft(self.k, fit_y5 - self.experiment * self.k ** self.weight )
        # plt.plot(r, amp, linewidth=self.linewidth, label='COSCd')

        r,amp,real,imag = calcfft(self.k, fit_y6 - self.experiment * self.k **3 )
        plt.plot(r, amp, linewidth=self.linewidth, label='COSCdCd')

        r,amp,real,imag = calcfft(self.k, fit_y7 - self.experiment * self.k **3 )
        plt.plot(r, amp, linewidth=self.linewidth, label='OOSCd')

        r,amp,real,imag = calcfft(self.k, fit_y8 - self.experiment * self.k **3 )
        plt.plot(r, amp, linewidth=self.linewidth, label='OOSCdCd')

        plt.xlim([0,5])
        plt.legend(fontsize=5)
        plt.title('FT of difference in nonlinear curve fit (scipy)')
        plt.xlabel('R')
        plt.ylabel('amp')
        plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/FT_fitNL_diff.pdf',format='pdf')
        plt.show()

        return fit_y1,fit_y2,fit_y3,fit_y6, fit_y7, fit_y8




    def LLE_model_OSCd(self, params):
        e = params[-1]
        exp = np.vstack((self.k_original,self.experiment_original)).transpose()
        k_new=np.sqrt(self.k**2-0.2625*e)
        # print(k_new)
        exp = self.interpolation(exp,k_new)
        model = self.EXAFS_model_OSCd(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(exp * self.k ** self.weight, loc=model))
        return LL

    def LLE_model_OSCd2(self, params):
        e = params[-1]
        exp = np.vstack((self.k_original,self.experiment_original)).transpose()
        k_new=np.sqrt(self.k**2-0.2625*e)
        # print(k_new)
        exp = self.interpolation(exp,k_new)
        model = self.EXAFS_model_OSCd2(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(exp * self.k ** self.weight, loc=model))
        return LL

    def LLE_model_O2S(self, params):
        e = params[-1]
        exp = np.vstack((self.k_original,self.experiment_original)).transpose()
        k_new=np.sqrt(self.k**2-0.2625*e)
        # print(k_new)
        exp = self.interpolation(exp,k_new)
        model = self.EXAFS_model_O2S(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(exp * self.k ** self.weight, loc=model))
        return LL
    #
    # def LLE_model_O2SCd2(self, params):
    #     model = self.EXAFS_model_O2SCd2(self.k, *params[:-1])
    #     LL = -np.sum(stats.norm.logpdf(self.experiment * self.k ** self.weight, loc=model))
    #     return LL

    def LLE_model_OS(self, params):
        e = params[-1]
        exp = np.vstack((self.k_original,self.experiment_original)).transpose()
        k_new=np.sqrt(self.k**2-0.2625*e)
        # print(k_new)
        exp = self.interpolation(exp,k_new)
        model = self.EXAFS_model_OS(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(exp * self.k ** self.weight, loc=model))
        return LL

    def LLE_model_SCd(self, params):
        e = params[-1]
        # print(*params)
        exp = np.vstack((self.k_original,self.experiment_original)).transpose()
        k_new=np.sqrt(self.k**2-0.2625*e)
        # print(k_new)
        # print('e',e)
        # print('knew',len(k_new),k_new)
        # print('exp',len(exp),exp)
        exp = self.interpolation(exp,k_new)
        model = self.EXAFS_model_SCd(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(exp * self.k ** self.weight, loc=model))
        return LL

    def LLE_model_COS(self, params):
        e = params[-1]
        exp = np.vstack((self.k_original,self.experiment_original)).transpose()
        k_new=np.sqrt(self.k**2-0.2625*e)
        exp = self.interpolation(exp,k_new)
        model = self.EXAFS_model_COS(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(exp * self.k ** self.weight, loc=model))
        return LL

    def LLE_model_COSCd(self, params):
        e = params[-1]
        exp = np.vstack((self.k_original,self.experiment_original)).transpose()
        k_new=np.sqrt(self.k**2-0.2625*e)
        exp = self.interpolation(exp,k_new)
        model = self.EXAFS_model_COSCd(self.k, *params[:-1])
        LL = -np.sum(stats.norm.logpdf(exp * self.k ** self.weight, loc=model))
        return LL



    def fit_shell(self,shell,fit_num,n):
        bounds_OSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['O'], self.bounds_dic['e0']))
        params_OSCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['O'] + self.params_dic['e0']
        bounds_OSCd2 = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['Cd2'], self.bounds_dic['O'], self.bounds_dic['e0']))
        params_OSCd2 = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['Cd2'] + self.params_dic['O'] + self.params_dic['e0']
        bounds_O2S = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['O2'], self.bounds_dic['e0']))
        params_O2S = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['O2'] + self.params_dic['e0']
        bounds_O2SCd2 = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['Cd'], self.bounds_dic['O'], self.bounds_dic['O'], self.bounds_dic['e0']))
        params_O2SCd2 = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['Cd2'] + self.params_dic['O'] + self.params_dic['O'] + self.params_dic['e0']
        bounds_OS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['e0']))
        params_OS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['e0']
        bounds_SCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['e0']))
        params_SCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['e0']
        bounds_COS = np.vstack((self.bounds_dic['S'], self.bounds_dic['O'], self.bounds_dic['C'], self.bounds_dic['e0']))
        params_COS = self.params_dic['S'] + self.params_dic['O'] + self.params_dic['C'] + self.params_dic['e0']
        bounds_COSCd = np.vstack((self.bounds_dic['S'], self.bounds_dic['Cd'], self.bounds_dic['O'], self.bounds_dic['C'], self.bounds_dic['e0']))
        params_COSCd = self.params_dic['S'] + self.params_dic['Cd'] + self.params_dic['O'] + self.params_dic['C'] + self.params_dic['e0']
        # print(params_SCd)
        exp = np.vstack((self.k_original,self.experiment_original)).transpose()


        plt.subplot(fit_num,2,1+n*2)
        # print('bound',len(bounds_SCd),bounds_SCd)
        model_shell = eval('self.LLE_model_'+shell)

        bound_shell = eval('bounds_'+shell)
        LLE_fitresult = differential_evolution(model_shell,bounds=bound_shell,strategy='randtobest1exp', popsize=150 ,maxiter=10000)
        print(shell,': ', LLE_fitresult)
        e = LLE_fitresult.x[-1]
        k_new=np.sqrt(self.k**2-0.2625*e)
        exp_new = self.interpolation(exp,k_new)
        fit_y = eval('self.EXAFS_model_'+shell+'(self.k, *LLE_fitresult.x[:-1])')
        plt.plot(self.k, fit_y, linewidth=self.linewidth, label='{:s} shell'.format(shell))
        plt.plot(self.k,exp_new*self.k**2,'k.', markersize=0.5, label=self.sample)
        plt.legend(fontsize=5)
        plt.xlabel('K')
        plt.yticks([])
        plt.ylabel('Chi')
        plt.tight_layout()
        plt.subplot(fit_num,2,2+n*2)
        r,amp,real,img = calcfft(self.k,fit_y*self.k)
        plt.plot(r,amp, linewidth= self.linewidth, label = shell)
        r,amp,real,img = calcfft(self.k, exp_new * self.k**3)
        plt.plot(r,amp,'k.',label = self.sample,markersize=0.5)
        plt.xlim([0,5])
        plt.xlabel('R')
        plt.yticks([])
        plt.legend(fontsize=5)
        plt.title('Fourier Transform')
        plt.tight_layout()
        # plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/{:s}fit_{:s}.pdf'.format(self.sample,shell),format='pdf')
        # plt.show()

    def LLE_fit(self):

        x1 = 'SCd'
        x2 = 'OS'
        x3 = 'OSCd'
        x4 = 'OSCd2'
        x5 = 'O2S'
        fits = (x1,x2,x3,x4,x5)
        fit_num = len(fits)
        width_lib = [6,10]
        fig_width = width_lib[(np.abs(width_lib - self.k_max)).argmin()]
        plt.figure(figsize=(fig_width,3*fit_num))
        for n,shell in enumerate(fits):
            self.fit_shell(shell,fit_num,n)
        plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/{:s}fit.pdf'.format(self.sample),format='pdf')
        plt.show()


        # plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/FT_fitLLE_diff.pdf',format='pdf')
        # plt.show()
        return












