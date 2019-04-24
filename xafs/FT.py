####  fourier transform of experiment
from xafs import *


def ftwindow( x, xmin=None, xmax=None, dx=1, dx2=None,
             window='hanning', _larch=None, **kws):
    """
    create a Fourier transform window array.

    Parameters:
    -------------
      x:        1-d array array to build window on.
      xmin:     starting x for FT Window
      xmax:     ending x for FT Window
      dx:       tapering parameter for FT Window
      dx2:      second tapering parameter for FT Window (=dx)
      window:   name of window type

    Returns:
    ----------
    1-d window array.

    Notes:
    -------
    Valid Window names:
        hanning              cosine-squared taper
        parzen               linear taper
        welch                quadratic taper
        gaussian             Gaussian (normal) function window
        sine                 sine function window
        kaiser               Kaiser-Bessel function-derived window

    """

    VALID_WINDOWS = ['han', 'fha', 'gau', 'kai', 'par', 'wel', 'sin', 'bes']

    if window is None:
        window = VALID_WINDOWS[0]
    nam = window.strip().lower()[:3]
    if nam not in VALID_WINDOWS:
        raise RuntimeError("invalid window name %s" % window)

    dx1 = dx
    if dx2 is None:  dx2 = dx1
    if xmin is None: xmin = min(x)
    if xmax is None: xmax = max(x)

    xstep = (x[-1] - x[0]) / (len(x)-1)
    xeps  = 1.e-4 * xstep
    x1 = max(min(x), xmin - dx1/2.0)
    x2 = xmin + dx1/2.0  + xeps
    x3 = xmax - dx2/2.0  - xeps
    x4 = min(max(x), xmax + dx2/2.0)

    if nam == 'fha':
        if dx1 < 0: dx1 = 0
        if dx2 > 1: dx2 = 1
        x2 = x1 + xeps + dx1*(xmax-xmin)/2.0
        x3 = x4 - xeps - dx2*(xmax-xmin)/2.0
    elif nam == 'gau':
        dx1 = max(dx1, xeps)

    def asint(val): return int((val+xeps)/xstep)
    i1, i2, i3, i4 = asint(x1), asint(x2), asint(x3), asint(x4)
    i1, i2 = max(0, i1), max(0, i2)
    i3, i4 = min(len(x)-1, i3), min(len(x)-1, i4)
    if i2 == i1: i1 = max(0, i2-1)
    if i4 == i3: i3 = max(i2, i4-1)
    x1, x2, x3, x4 = x[i1], x[i2], x[i3], x[i4]
    if x1 == x2: x2 = x2+xeps
    if x3 == x4: x4 = x4+xeps
    # initial window
    fwin =  np.zeros(len(x))
    if i3 > i2:
        fwin[i2:i3] = np.ones(i3-i2)

    # now finish making window
    if nam in ('han', 'fha'):
        fwin[i1:i2+1] = np.sin((np.pi/2)*(x[i1:i2+1]-x1) / (x2-x1))**2
        fwin[i3:i4+1] = np.cos((np.pi/2)*(x[i3:i4+1]-x3) / (x4-x3))**2
    elif nam == 'par':
        fwin[i1:i2+1] =     (x[i1:i2+1]-x1) / (x2-x1)
        fwin[i3:i4+1] = 1 - (x[i3:i4+1]-x3) / (x4-x3)
    elif nam == 'wel':
        fwin[i1:i2+1] = 1 - ((x[i1:i2+1]-x2) / (x2-x1))**2
        fwin[i3:i4+1] = 1 - ((x[i3:i4+1]-x3) / (x4-x3))**2

    return fwin

def FT_chi(k, chi,dx=3):
    timestep = (k[1]-k[0])
    time = np.zeros(2048, dtype='complex128')
    time[int(round(k[0]/timestep)):int(round(k.shape[0]+k[0]/timestep))] = k
    # print(time)
    signal = np.zeros(2048, dtype='complex128')
    win = ftwindow(k,dx=dx)
    signal[int(round(k[0]/timestep)):int(round(k.shape[-1]+k[0]/timestep))] = chi * win


    dr= np.pi/(time.shape[-1]*timestep)
    r= np.linspace(0,dr*signal.shape[-1]-dr,signal.shape[-1])

    # print('size of r ', r.shape[-1])
    amp = (timestep / np.sqrt(np.pi)) * np.fft.fft(signal)
    return r, abs(amp)
    # print('size of amp ', amp.shape[-1])
    # plt.plot(r, abs(amp))
    # plt.xlim([0,5])
    # plt.show()


def plot_fit_FT(sample, fit_method='LLE', dx=6):
    if fit_method == 'LLE':
        y1,y2,y3,y4,y5,y6,y7 = sample.LLE_fit()
    if fit_method == 'NL':
        y1,y2,y3,y4,y5,y6,y7 = sample.NLcurve_fit()

    for no,i in enumerate([y1,y2,y3,y4,y5,y6,y7]):
        label = ['SCd','OS','OSCd','COS','COSCd','OSCdCd','OOSCd']
        r,amp = FT_chi(sample.k,i,dx=dx)
        plt.plot(r,amp, linewidth= sample.linewidth, label = label[no])
    r,amp = FT_chi(sample.k, sample.experiment* sample.k**sample.weight,dx=dx)
    plt.plot(r,amp,'k.',label = sample.sample,markersize=0.5)
    plt.xlim([0,5])
    plt.legend()
    plt.title('Fourier Transform of fitting result')
    plt.savefig('/Users/Sophia/ownCloud/PhD/Statistic Analysis/figure/FT_fit{:s}.pdf'.format(fit_method), format='pdf')
    plt.show()
