import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from numba_stats import truncnorm, truncexpon, norm, expon
from scipy.stats import skewnorm
from scipy.stats import rv_continuous
from iminuit import Minuit
from iminuit.cost import ExtendedUnbinnedNLL,ExtendedBinnedNLL
from scipy import special
from datetime import datetime
from pathlib import Path

mpl.rc('xtick', labelsize=15) 
mpl.rc('ytick', labelsize=15) 
mpl.rcParams["font.size"] = 20
mpl.rcParams["lines.linewidth"] = 3
mpl.rcParams["mathtext.default"] = 'regular'
mpl.rcParams['lines.markersize'] = 1

# ### limit on LSC energy and limits on NaI coincidence energy
# x_limit = (100, 3000)
# NaI_limit = [[[0], [2, [1700, 2600]]], [[0], [4, [1300, 1900]]]]

startTime = datetime.now()

# xLimHigh = [200, 500]
# xLimMidL = [150,400]
# xLimMidH = [150,450]
xLimLeft = [900,3100]
xLimRight = [500,2700]

xLimPlotting = [0,3000]

# xLimgauss = [150,500]


# class truncErfc_gen(rv_continuous):
#     ''' Class for a truncated complementary error function
#     a and b are bounds
#     m is shift
#     s is scale
#     Shift and scale do not affect bounds
#     Normalization of pdf did not work for scipy 1.5.2, but does for scipy 1.10.1
#     '''

#     def _argcheck(self, a, b, m, s):
#         return (a < b) and (s > 0)

#     def _get_support(self, a, b, m, s):
#         return a, b

#     def _basIntegErfc(self, a, b):  # Integral between a and b, neither shifted nor scaled: int_ab erfc(x) dx
#         if True:
#             return (np.exp(-a * a) - np.exp(-b * b)) / np.sqrt(np.pi) + b * special.erfc(b) - a * special.erfc(
#                 a)  # This analytical method is faster than numerical integration
#         else:
#             return integrate.quad(special.erfc, a, b)[0]  # Slow normalization

#     def _integErfc(self, a, b, m,
#                    s):  # Integral of shifted and scaled function but with fixed bounds: int_ab erfc((x-m)/s) dx
#         return s * self._basIntegErfc((a - m) / s, (b - m) / s)

#     def _pdf(self, x, a, b, m, s):
#         return special.erfc((x - m) / s) / self._integErfc(a, b, m, s)  # Faster normalization

#     def _cdf(self, x, a, b, m, s):
#         return self._integErfc(a, x, m, s) / self._integErfc(a, b, m, s)

# truncErfc = truncErfc_gen(name='truncErfc', momtype=1)

def ReadInFileNaIChannels(file):
    
    totalChannels = []
    with open(file) as f: 
        next(f)
        i = 0 
        for line in f:
            lines = line.split('\n')[0].split(';')
            totalChannels.append(int(lines[1]))
            
            if i == 32: #Loops through 32 times to guarantee that each channel is seen at least twice. 
                break
            i +=1
    Channels = list(set(totalChannels))   
    
    data = []
    for c in range(len(Channels)):
        data.append([ [] for _ in range(3)])
    # channels = [2*n for n in range(numChannel)]
    # channels = Channels
    
    
    
    with open(file) as f:
        next(f)
        for line in f:
            lines = line.split('\n')[0].split(';')
            
            for c in range(len(Channels)):
                if int(lines[1]) == Channels[c]:
                    data[c][0].append(float(lines[2]))
                    data[c][1].append(float(lines[3]))
                    data[c][2].append(lines[5])
    return data, Channels

# def erf_gauss (x, x_limit0, x_limit1, n_erf, n_gauss, loc, scale, sigma, mu):

#     return (n_erf+n_gauss, n_erf * truncErfc.pdf(x, x_limit0, x_limit1, loc, scale) + n_gauss * truncnorm.pdf(x, x_limit0, x_limit1, loc=mu, scale=sigma))

def exp_gauss_skew(x,n_gauss,n_exp,tau,sigma,mu_gauss,mu_exp,skew, xLim1,xLim2):
    return n_gauss+n_exp, (n_gauss*skewnorm.pdf(x,a = skew, loc = mu_gauss, scale = sigma) + n_exp * truncexpon.pdf(x, xLim1, xLim2, loc = mu_exp, scale = tau))

def exp_gauss_skew_CDF(x,n_gauss,n_exp,tau,sigma,mu_gauss,mu_exp,skew, xLim1,xLim2):
    return (n_gauss*skewnorm.cdf(x,a = skew, loc = mu_gauss, scale = sigma) + n_exp * truncexpon.cdf(x, xLim1, xLim2, loc = mu_exp, scale = tau))

def exp_gauss_skew_CDF(x,n_gauss,n_exp,tau,sigma,mu_gauss,mu_exp,skew, xLim1,xLim2):
    return (n_gauss*skewnorm.cdf(x,a = skew, loc = mu_gauss, scale = sigma) + n_exp * truncexpon.cdf(x, xLim1, xLim2, loc = mu_exp, scale = tau))

def expFunc(x, n_exp, mu_exp, tau, xLim1,xLim2):
    return n_exp, n_exp * truncexpon.pdf(x, xLim1, xLim2, loc = mu_exp, scale = tau)

def gauss_skew(x,n_gauss,sigma,mu_gauss,skew):
    return n_gauss*skewnorm.pdf(x,a = skew, loc = mu_gauss, scale = sigma)

def expSkewGaussFit(data, ax, xLim):
    nBins = 100
    binWidth = 5
    bins, step = np.linspace(xLim[0], xLim[1], nBins+1, retstep=True)
    counts, mplbins, patches = ax[1].hist(data, bins = bins, density = False, log = False, histtype = 'step', lw = 3, label = 'Data')
    maxInd = np.where(counts == max(counts[30:])) #Determine the index of the largest bin in the histogram. 
    # plt.show()
    c = ExtendedBinnedNLL(counts,bins, exp_gauss_skew_CDF)
    #c = ExtendedBinnedNLL(truncHist,truncbins, exp_gauss_skew_CDF)
    #0: n_gauss
    #1: n_exp
    #2: tau
    #3: Sigma
    #4: Mu_gauss
    #5: mu_exp
    #6: skew
    initFit = [1.7877e4,8e1,25.5,5.33e2,2e3,9.5,1]
    m = Minuit(c ,n_gauss = initFit[0], n_exp = initFit[1],tau = initFit[2], sigma = initFit[3], mu_gauss = initFit[4], mu_exp = initFit[5], skew = initFit[6],xLim1 = xLim[0],xLim2 = xLim[1])
    m.limits["n_gauss", "n_exp", "tau", "sigma", "mu_gauss", "mu_exp"] = (0, None)
    m.fixed["xLim1","xLim2"] = True
    # m.simplex() #Another fitting procedure. Less accurate but faster. Seems to make my chi^2 go from 0.8 - 1.3
    m.migrad()
    m.hesse()  
    print(m)  
    print(f'chi^2 / ndof = {m.fval/m.ndof}')
    # print(f'chi^2 = {m.fval}') #prints the chi_^2
    # print(f'ndof = {m.ndof}') #prints the number of degrees of freedom 

    return m, initFit


def PlotFit(data,Bckdata,runTimes,saveFigPath,FitResults, FitErrors,xLim):
    fig, ax = plt.subplots(2,1, figsize = (15,10))
    

    fitValues2, initParams2= expSkewGaussFit(data,ax, xLim)

    xRange = np.linspace(xLim[0], xLim[1], 1000)
    nBins = 100
    bins, step = np.linspace(xLim[0], xLim[1], nBins+1, retstep=True)

    curveSkew = step * exp_gauss_skew(xRange,fitValues2.values[0],fitValues2.values[1],fitValues2.values[2],fitValues2.values[3],fitValues2.values[4],fitValues2.values[5], fitValues2.values[6], xLim[0],xLim[1])[1]

    curveExp = step * expFunc(xRange,fitValues2.values[1],fitValues2.values[5],fitValues2.values[2], xLim[0],xLim[1])[1]
    curveGauss = step * gauss_skew(xRange,fitValues2.values[0],fitValues2.values[3],fitValues2.values[4],fitValues2.values[6])

    mean = fitValues2.values[4]+ fitValues2.values[3]*(fitValues2.values[6]/(np.sqrt(1+fitValues2.values[6]**2)))*np.sqrt(2/np.pi)
    meanErr = np.sqrt(fitValues2.errors[4]**2+(fitValues2.values[6]/np.sqrt(1+fitValues2.values[6]**2)*np.sqrt(2/np.pi)*fitValues2.errors[3])**2+(fitValues2.values[3]*np.sqrt(2/np.pi)*(1/(np.sqrt(1+fitValues2.values[6]**2))-fitValues2.values[6]/(2*(1+fitValues2.values[6]**2))**(3/2))*fitValues2.errors[6])**2)
    print(f"Mean: {mean}")
    print(f"Mean Error: {meanErr}")
    
    for i in range(0,7):
        FitResults[i].append(fitValues2.values[i])
        FitErrors[i].append(fitValues2.errors[i])
    FitResults[7].append(mean)
    FitErrors[7].append(meanErr)
    # print(fitValues2.values[4])
    
    

    ##############################################################################################################################
    ###  Attempt to just plot the fit over top of the full histogram. The issue is that the histograms have different binnings ###
    ##############################################################################################################################

    xRangePlotting = np.linspace(xLimPlotting[0], xLimPlotting[1], 1000)
    nBins = 100
    bins, step = np.linspace(xLimPlotting[0], xLimPlotting[1], nBins+1, retstep=True)
    dataCounts, dataBins = np.histogram(data,bins)
    BckCounts,BckBins = np.histogram(Bckdata,bins)
    
    normalize = False
    
    
    ax[0].hist(dataBins[:-1],bins=dataBins,density = normalize,weights=dataCounts/runTimes[0],histtype='step',label = f"{saveFigPath.split('/')[-1].split('_')[2]} Integral Spectrum")
    # ax[0].hist(BckBins[:-1],bins=BckBins,density = normalize,weights=BckCounts/runTimes[1],histtype='step',label = f"Background Integral Spectrum")
    
    # ax[0].hist(data, bins = bins, density = False, log = False, histtype = 'step', lw = 3, label = f"{saveFigPath.split('/')[-1].split('_')[2]} Integral Spectrum")
    # ax[0].hist(Bckdata, bins = bins, density = False, log = False, histtype = 'step', lw = 3, label = f"Background Integral Spectrum")
    
    

    # ax[1].hist(data, bins = bins, density = False, log = False, histtype = 'step', lw = 3, label = 'Data')

    ax[1].plot(xRange,curveSkew,label = 'Fit', linewidth = 4, color = 'C1')
    ax[1].plot([],[],' ',label = f"Fit $\chi$^2 = {round(fitValues2.fval,2)}")
    ax[1].plot([],[],' ',label = f"Fit ndof = {fitValues2.ndof}")
    ax[1].plot(xRange, curveExp, linestyle = "dashed", color = 'black', alpha = 0.5)
    ax[1].plot(xRange, curveGauss, linestyle = "dashed", color = 'black', alpha = 0.5)

    ax[0].set_xlabel(r'Integral (ADU $\mu$s)')
    ax[1].set_xlabel(r'Integral (ADU $\mu$s)')

    ax[0].set_ylabel('Counts/bin/s')
    ax[1].set_ylabel('Counts/bin')
    
    
    #0: n_gauss
    #1: n_exp
    #2: tau
    #3: Sigma
    #4: Mu_gauss
    #5: mu_exp
    #6: skew
    # mean = fitValues2.values[4]+ fitValues2.values[3]*(fitValues2.values[6]/(np.sqrt(1+fitValues2.values[6]**2)))*np.sqrt(2/np.pi)
    ax[1].axvline(FitResults[7][-1], label = "Mean",color = 'green')
    ax[1].axvspan(FitResults[7][-1]-FitErrors[7][-1],FitResults[7][-1]+FitErrors[7][-1], color = 'green', alpha = 0.5)
    # ax[1].axvline(fitValues2.values[4], color = 'red', label = "mu")


    ax[1].legend(loc = 'best')
    ax[0].legend(loc = 'best')
    
    plt.savefig(saveFigPath, dpi = 400,bbox_inches='tight')
    # plt.show()
    # plt.close()


    return FitResults, FitErrors



filePath = '/home/nick/PhD/KDK+/Large_LSC_testing/Vertical_scatter_geometry_v4/No_collimation/2025_05_25/'
file = 'SDataR_Large_LSC_vessel_Cs137_triple_coinc_Vertical_scatter_v4_no_col.CSV'

bckfilePath = '/home/nick/PhD/KDK+/Large_LSC_testing/Vertical_scatter_geometry_v4/2025_03_04/'
bckfile = 'SDataR_Large_LSC_vessel_NaI_bck_triple_coinc_Vertical_scatter_v4_40CG_LSC.CSV'

data, channels = ReadInFileNaIChannels(f'{filePath}/{file}')
bckdata, bckchannels = ReadInFileNaIChannels(f'{bckfilePath}/{bckfile}')

fileName = file.split('.CSV')[0]
bckfileName = bckfile.split('/')[-1].split('.CSV')[0]

saveFilePath = f"{filePath}/{fileName}/figures"
Path(f"{saveFilePath}").mkdir(parents=True, exist_ok=True)

MaxTime = (data[0][0][-1]-data[0][0][0])/1e3
BckMaxTime = (bckdata[0][0][-1]-bckdata[0][0][0])/1e3
runtimes = [MaxTime,BckMaxTime]

fitResults = [[],[],[],[],[],[],[],[]]
fitErrors = [[],[],[],[],[],[],[],[]]


for i in range(2):
    integrals = np.asarray(data[i][1],dtype = 'float32')
    bckintegrals = np.asarray(bckdata[i][1],dtype = 'float32')
    if i == 0:
        limits = xLimLeft
    else:
        limits = xLimRight

    fitResults,fitErrors = PlotFit(integrals,bckintegrals,runtimes,f'{saveFilePath}/{fileName}channel_{i}_Peak_fit.png',FitResults=fitResults,FitErrors=fitErrors,xLim = limits)

   
f = open(f'{filePath}/{file}_fit.fit', "w+")
f.write("PMT, n_gauss, n_exp, tau, sigma, mu_gauss, mu_exp, skew, Mean\n")
for i in range(2):
    if i == 0:
        f.write(f"Left PMT")
    else:
        f.write("Right PMT")
    f.write(f",{fitResults[0][i]},{fitResults[1][i]},{fitResults[2][i]},{fitResults[3][i]},{fitResults[4][i]},{fitResults[5][i]},{fitResults[6][i]},{fitResults[7][i]}\n")
    
f.close()


f = open(f'{filePath}/{file}_err.err', "w+")
f.write("PMT, n_gauss, n_exp, tau, sigma, mu_gauss, mu_exp, skew, Mean\n")
for i in range(2):
    if i == 0:
        f.write(f"Left PMT")
    else:
        f.write("Right PMT")
    f.write(f",{fitErrors[0][i]},{fitErrors[1][i]},{fitErrors[2][i]},{fitErrors[3][i]},{fitErrors[4][i]},{fitErrors[5][i]},{fitErrors[6][i]},{fitErrors[7][i]}\n")
    
f.close()
    


