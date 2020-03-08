from scipy import stats
from EXAFS_analysis.xafs import interpolation
import numpy as np
from scipy.optimize import differential_evolution,Bounds
import multiprocess as mp
import os
from matplotlib import pyplot as plt
import seaborn as sns


class fit():
    def __init__(self,k,chi):
        self.sample_k = k
        self.sample_chi = chi
        self.name = 'unknown'
        return

    def set_fit_range(self,start, end):
        self.start = start
        self.end = end
        self.k_fit = np.linspace(start, end, int((end - start) * 40))
        print(' fit range is reset between ',self.start,'  and  ', self.end)
        return

    def set_sample_name(self,name):
        self.name = str(name)

    def set_shells(self,shells):
        self.shells = shells
        self.bounds = np.vstack(([0]*len(shells)*3+[-10],[12,5,0.5]*len(shells)+[10])).T
        for shell in self.shells:
            shell.exp_phase = interpolation(np.vstack((shell.k,shell.exp_phase)).T,self.k_fit)
            shell.exp_feff = interpolation(np.vstack((shell.k,shell.exp_feff)).T,self.k_fit)
            shell.emfp = interpolation(np.vstack((shell.k,shell.emfp)).T,self.k_fit)

    def model_shell(self,N,R,del_ss,shell):
        return N * self.k_fit ** 2 \
               * np.abs(shell.exp_feff) / (self.k_fit * R ** 2) \
               * np.sin(2 * self.k_fit * R +  shell.exp_phase) \
               * np.exp(-2 * R / shell.emfp) \
               * np.exp(-2 * del_ss * self.k_fit ** 2)

    def cal_chi(self,params):
        calculated_chi = 0
        for i, shell in enumerate(self.shells):
            calculated_chi += self.model_shell(*params[(i * 3):(i * 3 + 3)],shell)
        return calculated_chi

    def LLE_model(self,params):
        e = params[-1]
        k_new = np.sqrt(self.k_fit ** 2 - 0.2625 * e)
        data = np.vstack((self.sample_k,self.sample_chi)).transpose()
        while True:
            try:
                self.chi_fit = interpolation(data,k_new)
                break
            except:
                self.set_fit_range(self.start + 0.5, self.end - 0.5)
                k_new = np.sqrt(self.k_fit ** 2 - 0.2625 * e)

        calculated_chi = self.cal_chi(params)
        LL = -np.sum(stats.norm.logpdf(self.chi_fit * self.k_fit ** 2, loc=calculated_chi))
        return LL

    def fit(self,seed):

        LLE_fitresult = differential_evolution(self.LLE_model, bounds=self.bounds, strategy='currenttobest1bin',
                                               maxiter=10000,seed=seed)
        params = LLE_fitresult.x
        LL= self.LLE_model(params)
        params = np.append(params,LL)
        print('fitted parameters  with seed ' + str(seed) + '>>>>>>>>>>>\n')
        print(params)
        return params

    def fit_loop(self,loop_num=1000):
        parameter = {}
        #########################
        #
        # use parallel computing
        #
        #>>>>>>>>>>>>>>>>>>>>>>>>
        pool = mp.Pool(mp.cpu_count())
        print('CPU number :  ',mp.cpu_count())
        print('loop start ..............')
        for i in range(loop_num):
            print('\n Fit experiment in loop  >>>>>  ',i)
            seed = np.random.randint(500) * i + np.random.randint(500)
            parameter[i] = pool.apply_async(self.fit,args=(seed,)).get()
        pool.close()
        pool.join()
        print('loop finished')
        self.params = np.array([parameter[i] for i in range(loop_num)])


        np.savetxt(os.path.join(os.getcwd() + '/result/parameters_' + self.name + '.csv'), self.params, delimiter=',',
                   header='N,  R,  delss, '* len(self.shells) + 'e, LL')

        return _plot_(self.params,self.shells,self.name)

    def save_result(self):

        LL = self.params[:,-1]
        best_params = self.params[np.where(self.params == np.min(LL))[0][0]]
        print(' Best parameters : \n',best_params,np.min(LL))
        fit_result = self.cal_chi(best_params)
        e = best_params[-2]
        k_new = np.sqrt(self.k_fit ** 2 - 0.2625 * e)
        exp = np.vstack((self.sample_k,self.sample_chi)).T
        exp_e_shifted = interpolation(exp,k_new)

        data = np.vstack((self.k_fit,exp_e_shifted,fit_result)).T
        np.savetxt(os.path.join(os.getcwd() + '/result/fit_result_' + self.name +'.csv'),data,delimiter=',',
                   header = 'k, experiment with shifted e, fit result * k**2')

        plt.figure()
        plt.plot(self.k_fit,exp_e_shifted * self.k_fit**2,'k-', label=self.name + 'experiment with shifted e')
        plt.plot(self.k_fit,fit_result,label='fit result')
        plt.plot(self.k_fit,fit_result - exp_e_shifted * self.k_fit**2, label = 'fit difference')
        plt.legend(loc='lower right')
        plt.xlabel('k ($\AA^{-1}$)')
        plt.ylabel('$k^2\chi$')
        plt.tight_layout()
        plt.savefig(os.path.join(os.getcwd() + '/result/fit_plot_' + self.name + '.pdf'),format='pdf')
        plt.close()









class _plot_(object):
    def __init__(self,params,shells,name):
        self.params = params
        self.shells = shells
        self.name = name

    def pdf_hist(self):
        LL = self.params[:,-1]
        weight = (LL - np.max(LL))/(np.min(LL)-np.max(LL))
        weight = weight/np.sum(weight)

        plt.figure(figsize=(15,(len(self.shells) + 1) * 2.5))
        sns.set()


        for i,shell in enumerate(self.shells):
            N,R,delss = self.params[:,i*3:i*3+3].T

            plt.subplot(len(self.shells), 3, i * 3 + 2)
            R_mean = np.mean(R)
            range = np.max(np.abs(R - R_mean)) * 0.95
            range_index = np.where(np.logical_and(R >= R_mean - range, R <= R_mean + range))
            plt.hist(R[range_index], weights=weight[range_index], bins=50)
            plt.xlabel('shell ' + shell.linked_atoms + ' bond length ($\AA$)')
            plt.tight_layout()

            plt.subplot(len(self.shells),3,i*3+1)
            plt.hist(N[range_index],weights=weight[range_index],bins=50)
            plt.xlabel('shell '+shell.linked_atoms+' coordination number')
            plt.tight_layout()

            plt.subplot(len(self.shells),3,i*3+3)
            plt.hist(delss[range_index],weights=weight[range_index],bins=50)
            plt.xlabel('shell '+ shell.linked_atoms +' relative Debye-Waller factor compared with bulk')
            plt.tight_layout()

        print('\n Finish plotting fit analysis with weight, saving figure ..........')
        plt.savefig(os.path.join(os.getcwd() + '/result/fit_'+ self.name +'_histogram.pdf'),format='pdf')
        plt.close()





