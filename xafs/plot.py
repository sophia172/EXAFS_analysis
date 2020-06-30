from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib import cm
import pandas as pd
import os
cmap = cm.get_cmap('tab10')

class visulisation():
    def __init__(self,filelist):
        self.params = {}
        for file in filelist:
            file_path = os.path.join(os.getcwd(), 'result', 'parameters_' + file + '.dat')
            self.params[file] = np.genfromtxt(file_path)

    def pdf_fit_mean(self):
        return

    def pdf_hist(self):
        max_shell = (max([self.params[sample].shape[1] for sample in self.params.keys()]) - 1)/3
        plt.figure(figsize=(15, (max_shell + 1) * 2.5))
        sns.set()

        for i in range(max_shell):
            for m in [1,0,2]:

                plt.subplot(max_shell, 3, i * 3 + m+1)
                for j,sample in enumerate(self.params.keys()):
                    df = pd.DataFrame(sample.params[:,m*3:m*3+3].T,columns=['N','R','delss'])
                    LL = sample.params[:,-1]
                    weight = (LL - np.max(LL)) / (np.min(LL) - np.max(LL))
                    weight = weight / np.sum(weight)

                    print('\n Start plotting fit analysis with weight in ', sample.name, '..........')

                    color = list(cmap(j * 2 / 10))
                    color[-1] = 0.6
                    color = tuple(color)


                    if m ==1:

                        R_mean = np.mean(df['R'])
                        R_range = np.max(np.abs(df['R'] - R_mean)) * 0.95
                        range_index = np.where(np.logical_and(df['R']>= R_mean - R_range, df['R'] <= R_mean + R_range))[0]
                        plt.hist(df['R'][range_index],color=color, bins=50,linewidth=0.03)

                    elif m == 0:
                        plt.hist(df['N'][range_index], color=color, bins=50,linewidth=0.03)

                    elif m == 2:
                        plt.hist(df['delss'][range_index], color = color,bins=50,linewidth=0.03)



                if m == 1:
                    plt.xlabel('shell ' + str(i + 1) + ' - $R$ ($\AA$)')
                    plt.yticks([])

                    if i ==0:
                        plt.xlim([2.3,2.5])
                    elif i == 2:
                        plt.xlim([3.5,4.3])
                    plt.tight_layout()
                elif m == 0:
                    plt.xlabel('shell ' + str(i + 1) + ' - $N$')
                    plt.yticks([])
                    plt.tight_layout()
                elif m == 2:
                    plt.xlabel('shell ' + str(i + 1) + ' - $\Delta\sigma^2$')
                    plt.yticks([])
                    plt.tight_layout()
                    if i ==0:
                        plt.xlim([0,0.08])
                    elif i == 2:
                        plt.xlim([0,0.028])








        print(' Finish plotting fit analysis with weight, saving figure ..........')
        plt.savefig(result_path + 'fit_QDs_histogram1.pdf',format='pdf')
        plt.show()
        plt.close()

