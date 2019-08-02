
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
import seaborn as sns
from scipy.stats import norm


# import pylab, os, numpy, getpass, wx
# import sys
import multiprocess as mp
###########################MAIN#######################################
###########################   Change for different system#######################################

parameter={}
k_fit = np.linspace(2,16.5,500)                                                   ###
R1= 2.40                                                                          ### need change
R2 = 3                                                                            ###

###########################MAIN#######################################

###########################MAIN#######################################


def interpolation(f, x):
    g = interp1d(f[:, 0], f[:, 1])
    return g(x.tolist())



###############################################################
##
#
#  Calculate Feff and phase from experimental data
#
#
###################################################################

def extract_amp_phase_lambda(exp_base_path,lambda_path,result_path,N=1,R=1,path_name='unnamed'):
    k_temp = np.linspace(2, 17, 60000)
    FEFF = np.genfromtxt(lambda_path)
    lambda1 = FEFF[:, (0, 5)]
    lambda1 = interpolation(lambda1, k_fit)
    path_chi = np.genfromtxt(exp_base_path)[:,(0,1)]
    # k_new = np.sqrt(k_temp ** 2 + 0.2625 * 13.8)
    path_chi = interpolation(path_chi,k_temp)/k_temp**2
    # Find peaks for amp
    # peak_index = find_peaks(abs(path_chi), prominence=0.01, distance=10)[0]
    # k_peak = k_temp[peak_index]
    # chi_peak = path_chi[peak_index]
    # spline_knot = int(len(k_peak) / 5)
    # t = [k_peak[spline_knot], k_peak[spline_knot * 2], k_peak[spline_knot * 3], k_peak[spline_knot * 4]]
    # spl = LSQUnivariateSpline(k_peak, chi_peak,t)
    # amp = spl(k_fit)
    #
    amp = np.abs(hilbert(path_chi))
    amp = np.vstack((k_temp,amp)).transpose()
    amp = interpolation(amp,k_fit)
    phase = np.unwrap(np.angle(hilbert(path_chi)))+np.pi/2 - 2*k_temp*R
    phase = np.vstack((k_temp, phase)).transpose()
    phase = interpolation(phase,k_fit)
    Feff = amp/np.exp(-2*R/lambda1)/N*k_fit*R**2



    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(FEFF[:,0],FEFF[:,2],label='FEFF calculated Feff')
    plt.plot(k_fit,Feff,label='experiment Feff  including MSRD')
    plt.title('foundamental experiement Feff extraction')
    plt.legend()
    plt.tight_layout()

    plt.subplot(2, 1, 2)
    plt.plot(FEFF[:, 0], FEFF[:, 1]+FEFF[:, 3], label='FEFF calculated phase')
    plt.plot(k_fit, phase, label='experiment phase')
    plt.title('foundamental experiement phase extraction')
    plt.legend()
    plt.tight_layout()
    plt.savefig(result_path+'extracted_{:s}.pdf'.format(path_name),
                format='pdf')
    plt.close()
    data = np.vstack((k_fit,Feff,phase)).transpose()
    np.savetxt(result_path+'extracted_{:s}.dat'.format(path_name),data,
               header='k, Feff including MSRD, phase')
    return Feff,phase,lambda1


###############################################################
##
#
#  spepcify the FEFF and result path and shell
#
#
###################################################################
def FEFF_exp(f):
    global result_path, FEFF,shell
    FEFF, result_path, shell = f()

###############################################################
##
#
#  General case for model function
#  devide shell besed on different sample
#
###################################################################
def model_shell(x,N,R,del_ss):
    if R < R1            :
        path = 1
    elif R >=R2:
        path = 3
    elif R >=R1 and R <R2:
        path = 2
    return N * x ** 2 \
           * np.abs( FEFF['Feff'+str(path)]) / (x * R ** 2) \
           * np.sin(2 * x * R +  FEFF['phase'+str(path)]) \
           * np.exp(-2 * R /  FEFF['lambda'+str(path)]) \
           * np.exp(-2 * del_ss * x ** 2)


###############################################################
##
#
#  Define fitted experiment
#
#
###################################################################

def Ge_experiment(file):
    data = np.genfromtxt(file)[:,(0,1,2)]
    global k,chi
    k = data[1:,0]
    chi = data[1:,1]/k**3

def CdS_experiment(experiment):
    if experiment == 'pyspline_M311':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_M311_Aug18_pyspline.dat'
        data = np.genfromtxt(path)[:,(2,-1)]
        exp = data[:,(0,1)]
    elif experiment == 'pyspline_M322':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_M322_Aug18_pyspline.dat'
        data = np.genfromtxt(path)[:,(2,-1)]
        exp = data[:,(0,1)]
    elif experiment == 'pyspline_bulk':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_bulk_Aug18_pyspline.dat'
        data = np.genfromtxt(path)[:,(2,-1)]
        exp = data[:,(0,1)]
    elif experiment == 'athena_R':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_R_Nov17.chik'
        data = np.genfromtxt(path)[:,(0,1)]
        exp = data[:,(0,1)]
    elif experiment == 'athena_R_12sh':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_R_Nov17.chiq'
        data = np.genfromtxt(path)[:,(0,1)]
        data[1:, 1] = data[1:, 1] / data[1:, 0] ** 2
        exp = data[:,(0,1)]
    elif experiment == 'athena_bulk_12sh':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_bulk_Aug18_athena_12sh.chiq'
        data = np.genfromtxt(path)[:,(0,1)]
        data[1:,1]=data[1:,1]/data[1:,0]**2
        exp = data[:,(0,1)]
    elif experiment == 'athena_M311':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_XAFS_CdK_chi_2018.txt'
        data = np.genfromtxt(path)[:,(0,2,3,4)]
        exp = data[:,(0,1)]
    elif experiment == 'athena_M311_12sh':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_M311_2018.chiq'
        data = np.genfromtxt(path)[:,(0,1)]
        data[1:,1]=data[1:,1]/data[1:,0]**2
        exp = data[:,(0,1)]
    elif experiment == 'athena_bulk':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_XAFS_CdK_chi_2018.txt'
        data = np.genfromtxt(path)[:,(0,2,3,4)]
        exp = data[:,(0,2)]
    elif experiment == 'athena_bulk_ref':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_ref.chik'
        data = np.genfromtxt(path)[:,(0,1)]
        exp = data[:,(0,1)]
    elif experiment == 'athena_M322':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_XAFS_CdK_chi_2018.txt'
        data = np.genfromtxt(path)[:,(0,2,3,4)]
        exp = data[:,(0,3)]
    elif experiment == 'athena_M322_12sh':
        path = '/data/home/apw399/EXAFS_analysis/experiment/CdS_M322_2018.chiq'
        data = np.genfromtxt(path)[:,(0,1)]
        data[1:,1]=data[1:,1]/data[1:,0]**2
        exp = data[:,(0,1)]

    global k,chi
    k = exp[:,0]
    chi = exp[:,1]
    return


###############################################################
##
#
#  Define bounds for fit
#
#
###################################################################
bounds_dic = {'CdS': np.array([(0, 6), (R1, R2-0.0001), (-0.00001, 0.01)]),
              'CdO': np.array([(0, 12), (2.2, R1-0.0001), (0.0001, 0.05)]),
              'CdCd': np.array([(0, 15), (R2, 4.5), (-0.001, 0.03)]),
              'Ge1': np.array([(1, 6), (1, R1-0.0001), (0, 0.05)]),
              'Ge2': np.array([(1, 15), (R1, R2-0.0001), (0, 0.1)]),
              'Ge3': np.array([(1, 15), (R2, 5), (0, 0.1)]),
              'e': np.array([(-20, 10)]),
              'S0': np.array([(0.6, 1)])}

def Ge_bounds():
    global bounds
    bounds = np.vstack((bounds_dic['Ge1'],bounds_dic['e']))
    return

def CdOS_bounds():                                                                    ###
    global bounds                                                                     ### need change
    bounds = np.vstack((bounds_dic['CdS'],bounds_dic['CdS'],bounds_dic['CdCd'],bounds_dic['e']))
    return                                                                            ###

###############################################################
##
#
#  spepcify the FEFF calculation
#
#
###################################################################
def Ge_extract():
    FEFF_path1 = '/Users/sophia/ownCloud/PhD/Simulation/Ge/feff0001.dat'
    FEFF_path2 = '/Users/sophia/ownCloud/PhD/Simulation/Ge/feff0002.dat'
    FEFF_path3 = '/Users/sophia/ownCloud/PhD/Simulation/Ge/feff0003.dat'
    exp_path1 = '/Users/sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi1.chiq'
    exp_path2 = '/Users/sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi2.chiq'
    exp_path3 = '/Users/sophia/ownCloud/PhD/Experiment_Analysis/Ge/c_ge12_chi3.chiq'
    result_path = '/Users/sophia/ownCloud/PhD/Statistic Analysis/Ge/result/'


    shell_num = 1


    Ge_FEFF = {}
    Ge_FEFF['Feff1'], Ge_FEFF['phase1'], Ge_FEFF['lambda1'] = \
        extract_amp_phase_lambda(exp_path1,FEFF_path1, result_path, N=4,R=2.45,path_name='Ge-Ge1')
    Ge_FEFF['Feff2'], Ge_FEFF['phase2'], Ge_FEFF['lambda2'] = \
        extract_amp_phase_lambda(exp_path2, FEFF_path2, result_path, N=12, R=4, path_name='Ge-Ge2')
    Ge_FEFF['Feff3'], Ge_FEFF['phase3'], Ge_FEFF['lambda3'] = \
        extract_amp_phase_lambda(exp_path3, FEFF_path3, result_path, N=12, R=4.69, path_name='Ge-Ge3')
    return Ge_FEFF,result_path,shell_num


def CdSO_extract_cluster():
    FEFF_path2 = '/data/home/apw399/EXAFS_analysis/source/feff_CdS.dat'
    FEFF_path3 = '/data/home/apw399/EXAFS_analysis/source/feff_CdCd.dat'
    FEFF_path1 = '/data/home/apw399/EXAFS_analysis/source/feff_CdO.dat'
    exp_path2 = '/data/home/apw399/EXAFS_analysis/source/CdS_bulk_Aug18_athena_1sh.chiq'
    exp_path3 = '/data/home/apw399/EXAFS_analysis/source/CdS_bulk_Aug18_athena_2sh.chiq'
    exp_path1 = '/data/home/apw399/EXAFS_analysis/source/CdO_ref_1sh.chiq'
    ### need change
    result_path = '/data/home/apw399/EXAFS_analysis/result/2CdS_CdCd/'                   ###
    ### need change
    shell_num = 3                                                                      ###
    CdSO_FEFF = {}
    CdSO_FEFF['Feff2'], CdSO_FEFF['phase2'], CdSO_FEFF['lambda2'] = \
        extract_amp_phase_lambda(exp_path2,FEFF_path2, result_path, N=4,R=2.5118,path_name='Cd-S')
    CdSO_FEFF['Feff3'], CdSO_FEFF['phase3'], CdSO_FEFF['lambda3'] = \
        extract_amp_phase_lambda(exp_path3, FEFF_path3, result_path, N=12, R=4.1016, path_name='Cd-Cd')
    CdSO_FEFF['Feff1'], CdSO_FEFF['phase1'], CdSO_FEFF['lambda1'] = \
        extract_amp_phase_lambda(exp_path1, FEFF_path1, result_path, N=6, R=2.3475, path_name='Cd-O')
    return CdSO_FEFF,result_path,shell_num



###############################################################
##
#
#  Define fit model, used in fit()
#
#
###################################################################
def LLE_model(params):
    e = params[-1]
    k_new = np.sqrt(k_fit ** 2 - 0.2625 * e)
    data = np.vstack((k,chi)).transpose()
    chi_fit = interpolation(data,k_new)
    model_LLE = 0
    for i in range(shell):
        model_LLE += model_shell(k_fit, *params[(i * 3):(i * 3 + 3)])

    LL = -np.sum(stats.norm.logpdf(chi_fit * k_fit ** 2, loc=model_LLE))
    return LL



###############################################################
##
#
#  fit
#
#
###################################################################

def fit(seed):
    if seed % 2 != 0:
        LLE_fitresult = differential_evolution(LLE_model, bounds=bounds, strategy='currenttobest1bin', maxiter=1000,
                                               seed=seed)
    else:
        LLE_fitresult = differential_evolution(LLE_model, bounds=bounds, strategy='best1bin', maxiter=100000, seed=seed)

    params = LLE_fitresult.x
    model_LLE = 0
    for i in range(shell):
        model_LLE += model_shell(k_fit, *params[i*3:i*3+3])
    LL = LLE_model(params)
    e = params[-1]
    k_new = np.sqrt(k_fit ** 2 - 0.2625 * e)
    data = np.vstack((k, chi)).transpose()
    chi_fit = interpolation(data, k_new)
    params = np.append(params,LL)
    return chi_fit, model_LLE, params


###############################################################
##
#
#  put fit result in a dictionary
#  1. i is looped sequence
#  2. experiment is the file name which will be recorded
#
###################################################################
def gen_fit_result(experiment,i,seed=8754563):
    parameter[experiment + '_chi_fit' + str(i)],parameter[experiment+'_model_LLE'+str(i)], \
    parameter[experiment+'_params'+str(i)]  = fit(seed)

###############################################################
##
#
#  1. run fit in a looped times
#  2. store all results
#  3. find best fit
#
###################################################################

def loop_min_LL(experiment,loop_num):
    gen_fit_result(experiment, 0)
    parameter[experiment] = np.vstack((np.reshape(parameter[experiment + '_chi_fit' + str(0)], (-1, 1)),
                                       np.reshape(parameter[experiment + '_model_LLE' + str(0)], (-1, 1)),
                                       np.reshape(parameter[experiment + '_params' + str(0)], (-1, 1))))
    parameter[experiment + '_k'] = k
    parameter[experiment + '_chi'] = chi
    global results
    results = parameter[experiment]
    #########################
    #
    # use parallel computing
    #
    #>>>>>>>>>>>>>>>>>>>>>>>>
    pool = mp.Pool(mp.cpu_count())
    print('CPU number :  ',mp.cpu_count())

    def cal_fit(experiment, i):
        print('################################ loop ', i)
        gen_fit_result(experiment, i, seed=i + np.random.randint(50))
        parameter[experiment + str(i)] = np.vstack((np.reshape(parameter[experiment + '_chi_fit' + str(i)], (-1, 1)),
                                                    np.reshape(parameter[experiment + '_model_LLE' + str(i)], (-1, 1)),
                                                    np.reshape(parameter[experiment + '_params' + str(i)], (-1, 1))))
        return parameter[experiment + str(i)]

    def collected_result(result):
        global results
        results = np.append(results, result, axis=1)

    for i in range(loop_num):
        parameter[experiment+str(i)] = pool.apply_async(cal_fit, args=(experiment, i),callback=collected_result)

    pool.close()
    pool.join()

    print('Finished Looping, Now saving result to data file.........')
    #<<<<<<<<<<<<<<<<<<<<<<
    #
    # use parallel computing
    #
    #############################
    parameter[experiment] = results
    min_LL_index = parameter[experiment][-1,:].argmin()

    parameter[experiment + '_chi_fit'] = parameter[experiment][:len(k_fit),min_LL_index]
    parameter[experiment + '_model_LLE'] =  parameter[experiment][len(k_fit):len(k_fit)*2,min_LL_index]
    parameter[experiment + '_LL'] = parameter[experiment][-1,min_LL_index]
    parameter[experiment + 'e'] = parameter[experiment][-2, min_LL_index]
    parameter[experiment + '_params_list'] = parameter[experiment][-shell*3-2:,:]
    parameter[experiment + '_params'] = parameter[experiment][-shell*3-2:-2,min_LL_index]

    np.savetxt(result_path+experiment+'_params_all.dat',parameter[experiment][-shell*3-2:,:].transpose(),header='N R delss e LL')

    return



###############################################################
##
#
#  find the best fit in all fitted results and plot
#
#
###################################################################

def plot_best_fit(experiment, loop_num=1000):
    loop_min_LL(experiment, loop_num)
    print('\n Start plotting best fit in ', experiment, '..........')
    plt.figure()

    plt.plot(parameter[experiment + '_k'], parameter[experiment + '_chi'] * parameter[experiment + '_k'] ** 2, '.', markersize=0.8,
             label=experiment)
    plt.plot(k_fit, parameter[experiment + '_chi_fit'] * k_fit ** 2, linewidth=0.5, label=experiment+'_e_shifted')
    plt.plot(k_fit, parameter[experiment + '_model_LLE'], linewidth=0.5, label=experiment + '_fit')
    plt.plot(k_fit, parameter[experiment + '_model_LLE'] - parameter[experiment + '_chi_fit'] * k_fit ** 2,
             linewidth=0.5, label=experiment + '_difference')
    data = np.vstack((k_fit, parameter[experiment + '_chi_fit'],
                      parameter[experiment + '_model_LLE'] / k_fit ** 2,
                      parameter[experiment + '_model_LLE'] / k_fit ** 2
                      - parameter[experiment + '_chi_fit'])).transpose()
    np.savetxt(result_path + experiment + '_fit.dat', data, header='k    chi_fitted with shifted e      fit/chi     difference')
    plt.legend(loc='lower right')
    plt.xlabel('k ($\AA^{-1}$)')
    plt.ylabel('$k^2\chi$')
    plt.tight_layout()
    print(parameter[experiment + '_params'])
    print(shell)
    params_text = ['shell'+ str(i+1) +'  {:.2f}  {:.3f}  {:.5f} \n'.format(*parameter[experiment + '_params'][3*i:3*i+3])
                   for i in range(shell) ]
    params_label = ''
    for i in range(shell):
        params_label += params_text[i]

    plt.text(9, 0.2, 'LL = ' + str(parameter[experiment + '_LL']) + '\n'
             + ' e  = '  + str(parameter[experiment + 'e']) + '\n'
             + '            N        R      del_ss\n'
             + params_label)

    print(' Finish plotting best fit, Save figure ........')
    plt.savefig(result_path + 'fit_' + experiment + '_(extracted_feature).pdf',format='pdf')
    plt.close()
    return

###############################################################
##
#
#  statistic distribution of fit result
#
#
###################################################################

def plot_fit_hist(experiment):

    def param_list(i):
        N = parameter[experiment+'_params_list'][i*3,:]
        R = parameter[experiment+'_params_list'][i*3+1,:]
        delss = parameter[experiment+'_params_list'][i*3+2,:]
        return N,R,delss

    LL = parameter[experiment+'_params_list'][-1,:]
    weight = (LL - np.max(LL))/(np.min(LL)-np.max(LL))
    weight = weight/np.sum(weight)

    print('\n Start plotting fit analysis with weight in ', experiment, '..........')
    plt.figure(figsize=(15,(shell+1)*2.5))
    sns.set()

    for i in range(shell):
        N,R,delss = param_list(i)
        plt.subplot(shell,3,i*3+1)
        plt.hist(N,weights=weight,bins=50)
        plt.xlabel('shell '+str(i+1)+' coordination number')
        plt.tight_layout()

        plt.subplot(shell,3,i*3+2)
        plt.hist(R,weights=weight,bins=50)
        plt.xlabel('shell '+str(i+1)+' bond length ($\AA$)')
        plt.tight_layout()

        plt.subplot(shell,3,i*3+3)
        plt.hist(delss,weights=weight,bins=50)
        plt.xlabel('shell '+str(i+1)+' relative Debye-Waller factor compared with bul')
        plt.tight_layout()

    print(' Finish plotting fit analysis with weight, saving figure ..........')
    plt.savefig(result_path + 'fit_'+experiment+'_histogram.pdf',format='pdf')
    plt.close()

    weight = np.ones(np.shape(LL)[0])/np.shape(LL)[0]

    print('\n Start plotting fit analysis in ', experiment, ' without weight..........')
    plt.figure(figsize=(15,(shell+1)*2.5))
    sns.set()

    for i in range(shell):
        N,R,delss = param_list(i)
        plt.subplot(shell,3,i*3+1)
        plt.hist(N,weights=weight,bins=50)
        plt.xlabel('shell '+str(i+1)+' coordination number')
        plt.tight_layout()

        plt.subplot(shell,3,i*3+2)
        plt.hist(R,weights=weight,bins=50)
        plt.xlabel('shell '+str(i+1)+' bond length ($\AA$)')
        plt.tight_layout()

        plt.subplot(shell,3,i*3+3)
        plt.hist(delss,weights=weight,bins=50)
        plt.xlabel('shell '+str(i+1)+' relative Debye-Waller factor compared with bul')
        plt.tight_layout()

    print('Finish plotting fit analysis without weight, saving figure ..........')
    plt.savefig(result_path + 'fit_'+experiment+'_histogram(no_weight).pdf',format='pdf')
    plt.close()
    return

#############################################################
#
# calculate covariance and correlation matrix here
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def correlation(experiment):
    print('\n Analyzing correlation in fitting parameter in ', experiment, '......... ')
    variable = parameter[experiment][-3*shell-2:-2,:]
    covariance = np.cov(variable)
    corr,_=stats.spearmanr(variable,axis=1)
    return corr
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
# calculate covariance and correlation matrix here
#
# ############################################################



#############################################################
#
# plot correlation
#
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
def plot_corr(experiment):
    corr = correlation(experiment)
    sns.set(style="white")
    mask = np.zeros_like(corr, dtype=np.bool)
    mask[np.triu_indices_from(mask)] = True
    f, ax = plt.subplots(figsize=(11, 9))
    cmap = sns.diverging_palette(220,10, as_cmap=True)
    label = []
    for i in range(shell):
        label.append('shell'+ str(i+1) +'_N')
        label.append('shell'+ str(i+1) +'_R')
        label.append('shell'+ str(i+1) +'_delss')

    print('\n Start plotting correlation analysis in ', experiment, '..........')
    sns.heatmap(corr, mask=mask, cmap=cmap, vmax=1, annot=True, center=0,
                square=True, linewidths=.5, cbar_kws={"shrink": .5},
                xticklabels=label,
                yticklabels=label)
    print('\n Finish plotting correlation analysis, saving figure..........')
    plt.savefig(
        result_path + '/correlation_' + experiment + '.pdf',
        format='pdf')
    plt.close()




############################################################
#
# Fit for Gemanium experiment
#
############################################################
# filelist = glob.glob('/Users/sophia/ownCloud/PhD/Experiment_Analysis/Ge/Ge_QDs_HP_Nano_Letters/Ge_HP_Chi-k/*.k3')
#     FEFF_exp(Ge_extract)
# for file in filelist:
#     filename = os.path.basename(file)
#     Ge_experiment(file)
#     Ge_bounds()
#     plot_best_fit(filename[:-3], loop_num=10)
#     plot_fit_hist(filename[:-3])
#     plot_corr(filename[:-3])
#     print(filename)


############################################################
#
# Fit for CdS experiment
#
############################################################
experiment_list = ['athena_bulk_12sh','athena_bulk','athena_R','athena_M311','athena_M322','athena_R_12sh','athena_M311_12sh','athena_M322_12sh','pyspline_M311','pyspline_M322','pyspline_bulk']
FEFF_exp(CdSO_extract_cluster)
CdOS_bounds()
def parel_plot(experiment,loop_num=10):
    print('\n \n \n Start processing experiment data :  ', experiment,' .......................')
    CdS_experiment(experiment)
    plot_best_fit(experiment,loop_num=loop_num)
    plot_fit_hist(experiment)
    plot_corr(experiment)
    return
for i in experiment_list[:6]:
    parel_plot(i,loop_num=10000)

