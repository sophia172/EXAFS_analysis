####  fourier transform of experiment
from xafs import *

DR=0.05 #step size of R
fftPOINTS=1048 #number of fft points; works best if equal to 2^n




def calcfft(kdata,k3xafsdata,kmin=3,kmax=19,windowType='hanning'):
    #calculate dk
    dk=pi/(fftPOINTS*DR)
    # 1/sqrt(pi)
    sqrtpi = 0.56418958350
    # Normalization
    cnorm = dk * sqrtpi
    index=0
    # Bin the data to a grid
    bindata=[]
    krange=arange(dk,kdata[-1],dk)
    kdata = kdata.tolist()
    for bin in krange:
        mysum=0
        for k in kdata[index:]: #only count from last bin
            if(k<bin): #if k falls into bin
                last=kdata.index(k)
                temp=k3xafsdata[last]
                mysum+=temp
            else:
                last=kdata.index(k) #first k that doesn't fall in bin
                break
        bindata.append(mysum*cnorm/(last-index+1))
        index=last

    # apply window
    index_kmin = (np.abs(krange - kmin)).argmin()
    index_kmax = (np.abs(krange - kmax)).argmin()


    window_len=index_kmax-index_kmin+1
    if not windowType in ['flat', 'hanning', 'hamming', 'bartlett','blackman']:
        windowType=='flat'
    if windowType == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+windowType+'(window_len)')
    for i in range(len(bindata)):

        if(i>=index_kmin and i<=index_kmax):
            bindata[i]=bindata[i]*w[i-index_kmin]
        else:
            bindata[i]=0.0


    # FFT
    rawfftdata=fft.fft(bindata,fftPOINTS)
    conjfftdata=conjugate(rawfftdata)
    fftdata=sqrt(rawfftdata*conjfftdata).real
    fftlist=fftdata.tolist()
    data=[]
    realdata=[]
    imagdata=[]
    for i in range(int(len(fftlist)/2)):
        data.append((fftlist[i]+fftlist[-i])/2.0)
        realdata.append((rawfftdata.real[i]+rawfftdata.real[-i])/2.0)
        imagdata.append((rawfftdata.imag[i]-rawfftdata.imag[-i])/2.0)

    # Return magnitude, real and imaginary data
    r = arange(0,len(data)*DR,DR)
    return r,data,realdata,imagdata


def calcifft(chir, nfft=2048):
    """
    calculate reverse XAFS Fourier transform, from chi(R) to
    chi(q), using common XAFS conventions.  This version demands
    chir be the complex chi(R) as created from xftf().

    It returns the complex array of chi(q) without putting any
    values into an output group.

    Parameters:
    -------------
      chir:     1-d array of chi(R) to be transformed
      nfft:     value to use for N_fft (2048).
      kstep:    value to use for delta_k (0.05).

    Returns:
    ----------
      complex 1-d array for chi(q).

    This is useful for repeated FTs, as inside loops.
    """
    kstep= np.pi/(fftPOINTS*DR)
    cchi = zeros(nfft, dtype='complex128')
    cchi[0:len(chir)] = chir
    return  (4*sqrt(pi)/kstep) * np.fft.ifft(cchi)[:int(nfft/2)]

def plot_fit_FT(sample, fit_method='LLE'):
    if fit_method == 'LLE':
        y1,y2,y3,y4,y5 = sample.LLE_fit()
    if fit_method == 'NL':
        y1,y2,y3,y4,y5,y6 = sample.NLcurve_fit()
    plt.subplot(2,1,2)
    for no,i in enumerate([y1,y2,y3,y4,y5]):
        label = ['SCd','OS','OSCd','OSCdCd','OOS']
        r,amp,real,img = calcfft(sample.k,i*sample.k)
        plt.plot(r,amp, linewidth= sample.linewidth, label = label[no])
    r,amp,real,img = calcfft(sample.k, sample.experiment * sample.k**3)
    plt.plot(r,amp,'k.',label = sample.sample,markersize=0.5)
    plt.xlim([0,5])
    plt.xlabel('R')
    plt.yticks([])
    plt.legend(fontsize=5)
    # plt.title('Fourier Transform of fitting result')
    plt.tight_layout()
    plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/fit{:s}.pdf'.format(fit_method), format='pdf')
    plt.show()
