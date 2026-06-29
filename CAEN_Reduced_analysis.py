import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from pathlib import Path
import scipy.stats as Stats
from iminuit import Minuit
from iminuit.cost import ExtendedUnbinnedNLL, ExtendedBinnedNLL
from numba_stats import truncnorm, truncexpon
from scipy.stats import skewnorm
import matplotlib

matplotlib.rcParams["font.size"] = 20
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["mathtext.default"] = 'regular'
matplotlib.rcParams['lines.markersize'] = 3


class ChData:
    def __init__(self,ch,E,t,flag):
        self.ch = ch
        self.E = [E]
        self.t = [t]
        self.flag = [flag]
        
        self.fitWindow = None
        
    def AddEvent(self,E,t,flag):
        self.E.append(E)
        self.t.append(t)
        self.flag.append(flag)
        
    def EnergyHist1DPlot(self, ax, detectors,ChBins, BinRange, norm = False, log = False):
        #############################################################################
        #   This is a function that plots the energy spectra of a given event.      #
        #   The inputs of this function are the detector object itself, the axis    #
        #   that it will be plotted over, then the number of bins and the bin       #
        #   range. By default the spectra are all normalized to rate. I have no     #
        #   plans to change this right now, but maybe in the future it will be      #
        #   changed to be an option.                                                #
        #############################################################################
        #   Variables:                                                              #
        #       - self: The detector object that will be plotted.                   #
        #       - ax: The axis from plt.subplots that this spectra will be plotted  #
        #             on.                                                           #
        #       - ChBins (int): The number of bins that will be used for this       #
        #                       histogram.                                          #
        #       - BinRange (array(int)): The range that the data will be plotted    #
        #                                over.                                      #
        #       - norm (boolean): Normalizes the data to an area of 1.              #
        #       - log (boolean): Puts the y-axis to log scale.                      #
        #############################################################################
        
        self.Ehist, self.EhistBins = BinHistograms(self.E,BinRange,ChBins) #Calls the bin histogram function to bin the data. 
        ax.hist(self.EhistBins[:-1],bins =self.EhistBins,density = norm,weights=self.Ehist/(self.t[-1] - self.t[0]),histtype = 'step',label = f'{detectors[f'Channel {self.ch}'][0]}', log = log,linewidth = 3) #Plot the data on the provided histogram. Note that this data is always rate normalized. 
            
        
        ax.plot([],[],' ', label = f'Events: {len(self.E)}')
        ax.plot([],[],' ', label = f'Run Time: {round((self.t[-1]-self.t[0])/1e9,3)} s')
        #Set various labels for the plots. 
        ax.set_title(detectors[f'Channel {self.ch}'][0])
        ax.set_xlabel('Integral (ADC Channel)')
        ax.set_ylabel('Counts/bin/ns')
        ax.legend(loc = 'best')
        
    def EnergyStabilityPlot(self, ax,fig,detectors, TimeBinRange,ChBinRange):
        
        h = ax.hist2d(self.t, self.E, bins = (TimeBinRange,ChBinRange),cmin = 1)
        fig.colorbar(h[3],ax=ax)
        
        ax.set_title(detectors[f'Channel {self.ch}'][0])
        ax.set_xlabel('Trigger Time (ns)')
        ax.set_ylabel('Integral (ADC) Channel')
        
    def FlagStabilityPlot(self,ax,detectors):
        
        uniqueFlags = sorted(set(self.flag))
        
        flagDictionary = dict(zip(uniqueFlags, [0 for _ in range(len(self.flag))]))
        
        for i in (self.flag):
            flagDictionary[i] +=1
            
        ax.bar(range(len(flagDictionary)), list(flagDictionary.values()),align = 'center')
        ax.set_title(f'{detectors[f'Channel {self.ch}'][0]} flag distribution')
        ax.set_yscale('log')
        ax.set_xticks(range(len(flagDictionary)), list(flagDictionary.keys()),rotation = 70)
        ax.set_ylabel('Number of Events')
        
    def FlagHistograms(self,ax,detectors):
        
        uniqueFlags = sorted(set(self.flag))
        
        flagDictionary = dict(zip(uniqueFlags, [[] for _ in range(len(self.flag))]))

        for i,flag in enumerate(self.flag):
            flagDictionary[flag].append(self.E[i])
        
        for flag in flagDictionary:
            flagDictionary[flag]
            
            Ehist, EhistBins = BinHistograms(flagDictionary[flag],[0,4050],100) #Calls the bin histogram function to bin the data. 
            ax.hist(EhistBins[:-1],bins =EhistBins,weights=Ehist/(self.t[-1] - self.t[0]),histtype = 'step',label = flag,linewidth = 3)
            #(ax = ax[0], ChBins =ChBins[i], BinRange = Binrange, norm = False, log = False) 
        ax.set_title(f'{detectors[f'Channel {self.ch}'][0]} flag histograms')
        ax.set_yscale('log')
        ax.set_xlabel('Integral (ADU)')
        ax.set_ylabel('counts/bin/ns')
        # Shrink current axis by 20%
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        
    #####################################################################################################################################
    #  I've put a fit function here, but I am not sure that I will be using this that much. I will leave it for now.                    #
    #####################################################################################################################################  
        
    def expGaussFit(self, ax, xLim):
        #########################################################################
        #   Function that will fit the data that has been provided, as well as  #
        #   plotting the data for the histogram that was fit.                   #
        #########################################################################
        #   Variables:                                                          #
        #       - data: unbinned histogram data.                                #
        #       - ax: matplotlib axis that I will plot the histogram on.        #
        #       - xLim: The range over which this data will be plotted.         #
        #########################################################################
        nBins = 100
        binWidth = 5
        bins, step = np.linspace(xLim[0], xLim[1], nBins+1, retstep=True)
        counts, mplbins, patches = ax.hist(self.E, bins = bins, density = False, log = False, histtype = 'step', lw = 3, label = 'Data')
        # plt.show()
        maxInd = np.where(counts == max(counts[30:])) #Determine the index of the largest bin in the histogram. 
        c = ExtendedBinnedNLL(counts,bins, exp_gauss_CDF)
        #c = ExtendedBinnedNLL(truncHist,truncbins, exp_gauss_skew_CDF)
        #0: n_gauss
        #1: n_exp
        #2: tau
        #3: Sigma
        #4: Mu_gauss
        #5: mu_exp
        self.n_gauss = sum(counts)
        self.n_exp = 0.8*self.n_gauss
        self.tau = 100
        self.Mu_gauss = (xLim[1]+xLim[0])/2 #counts[maxInd[0]]
        self.Mu_exp = xLim[0]
        self.sigma = 0.1*self.Mu_gauss
        
        initFit = [self.n_gauss,self.n_exp,self.tau,self.sigma,self.Mu_gauss,self.Mu_exp] #NaI raw fit params
        # initFit = [9e3,100,50,200,1800,200,1]
        
        m = Minuit(c ,n_gauss = initFit[0], n_exp = initFit[1],tau = initFit[2], sigma = initFit[3], mu_gauss = initFit[4], mu_exp = initFit[5],xLim1 = xLim[0],xLim2 = xLim[1])
        m.limits["n_gauss", "n_exp", "tau", "sigma", "mu_gauss", "mu_exp"] = (0, None)
        m.fixed["xLim1","xLim2"] = True
        # print(m)
        #m.simplex() #Another fitting procedure. Less accurate but faster. Seems to make my chi^2 go from 0.8 - 1.3
        m.migrad()
        m.hesse()    
        print(r'$chi$^2 / ndof = ' + f"{m.fval/m.ndof}")
        print(m)
        # print(m.values[4])
        # print(f'chi^2 = {m.fval}') #prints the chi_^2
        # print(f'ndof = {m.ndof}') #prints the number of degrees of freedom 
        
        # self.mean = m.values[4] + m.values[3] * (m.values[6]/(np.sqrt(1+m.values[6]**2)))*np.sqrt(2/np.pi)
        self.mean = m.values[4]
        
        xRange = np.linspace(xLim[0], xLim[1], 1000)
        expGauss = step * exp_gauss(xRange, m.values[0],m.values[1],m.values[2],m.values[3],m.values[4],m.values[5], xLim[0], xLim[1])[1]
        exp = step * expFunc(xRange,m.values[1],m.values[5],m.values[2], xLim[0], xLim[1])[1]
        gaussFit = step * gauss(xRange,m.values[0], m.values[3],m.values[4])
        
        # initexpGauss = step * exp_gauss_skew(xRange, initFit[0],initFit[1],initFit[2],initFit[3],initFit[4],initFit[5],initFit[6], xLim[0], xLim[1])[1]
        # initexp = step * expFunc(xRange,initFit[1],initFit[5],initFit[2], xLim[0], xLim[1])[1]
        # initgauss = step * gauss_skew(xRange,initFit[0], initFit[3],initFit[4],initFit[6])
        
        ax.plot(xRange,expGauss, label = 'fit', linewidth = 4)
        ax.plot(xRange,exp, linestyle = "dashed", color = 'black', alpha = 0.5)
        ax.plot(xRange,gaussFit, linestyle = "dashed", color = 'black', alpha = 0.5)
        ax.axvline(self.mean,color = 'red', label = f'Mean = {round(self.mean,2)} +/- {round(m.errors[4],2)}')
        ax.plot([],[],' ', label = r"$\chi$^2 = " + f"{round(m.fval,2)}")
        ax.plot([],[],' ', label = r"dof = " + f"{round(m.ndof,2)}")
        ax.plot([],[],' ', label = r"$\chi$^2/dof = " + f"{round(m.fval/m.ndof,2)}")
        ax.legend(loc = 'best')
        # ax.plot([],[],' ',label = f"Fit ndof = {m.ndof}")
        
        
        return m
        
        
class coincData:
    def __init__(self,channel,time,Energy,flag):
        self.chList = channel
        self.chData = {}
        for i,ch in enumerate(channel):
            self.chData.update({ch : ChData(int(ch),Energy[i],time[i],flag[i])})
        
        # print(self.chData['4'].ch)
        
    def AddEvent(self,channel,time,Energy,flag):
        for i, ch in enumerate(channel):
            self.chData[ch].AddEvent(Energy[i],time[i],flag[i])
            
            
    def Stabilityplots(self,ChBins,saveFilePath,fileName,scale,detectors):
        
        
        
        fig, ax = plt.subplots(4,len(self.chList), figsize = (10*len(self.chList),25))
        plt.tight_layout()
        plt.subplots_adjust(left = 0.06,wspace = 0.2,hspace = 0.35,top = 0.98,bottom = 0.04)
        
        for i,coinc in enumerate(self.chData):
            Binrange = [0,4050]
            
            ChBinRange = np.linspace(Binrange[0],Binrange[1],ChBins[i])
            
            timebinsRange = np.linspace(min(self.chData[coinc].t),max(self.chData[coinc].t),100)
            
            if len(self.chList) ==1:
                self.chData[coinc].EnergyHist1DPlot(ax = ax[0], ChBins =ChBins[i], BinRange = Binrange, norm = False, log = False,detectors = detectors) 
            
                self.chData[coinc].EnergyStabilityPlot(ax[1],fig,TimeBinRange = timebinsRange, ChBinRange = ChBinRange,detectors = detectors)
                
                self.chData[coinc].FlagStabilityPlot(ax[2],detectors = detectors)
                self.chData[coinc].FlagHistograms(ax[3],detectors = detectors)
            else:
                self.chData[coinc].EnergyHist1DPlot(ax = ax[0,i], ChBins =ChBins[i], BinRange = Binrange, norm = False, log = False,detectors = detectors) 
                
                self.chData[coinc].EnergyStabilityPlot(ax[1,i],fig,TimeBinRange = timebinsRange, ChBinRange = ChBinRange,detectors = detectors)
                
                self.chData[coinc].FlagStabilityPlot(ax[2,i],detectors = detectors)
                self.chData[coinc].FlagHistograms(ax[3,i],detectors = detectors)
                
            
        if scale:
            plt.savefig(saveFilePath/f'{fileName}_scaled_stability_plots.png',bbox_inches='tight')
        else:
            plt.savefig(saveFilePath/f'{fileName}_stability_plots.png',bbox_inches='tight')
            
    def fitData(self,ChBins, Binrange, saveFilePath, fileName,fitWindow, detectors):
        
        
            
        # ChBinRange = np.linspace(Binrange[0],Binrange[1],ChBins[i])
        fig,ax = plt.subplots(2,len(self.chList), figsize = (15*len(self.chList),20))
        
        plt.tight_layout()
        plt.subplots_adjust(hspace = 0.3, wspace = 0.15, left = 0.1, right = 0.95, top = 0.95, bottom = 0.07)
        
        f = open(saveFilePath/f'{fileName}_fit_results.txt',"w+")
        f.write("Detector, n_gauss, n_exp, tau, sigma, mu_gauss, mu_exp, skew,$chi$^2 / ndof \n")
        for i,data in enumerate(self.chData): #Loops through all the detector objects. Then plot the 1D hist and the fit.
            self.chData[data].fitwindow = fitWindow[i]

            if len(self.chList) == 1:
                self.chData[data].EnergyHist1DPlot(ax = ax[0], ChBins =ChBins, BinRange = Binrange[i], norm = False, log = False,detectors = detectors)#Plots 
                m = self.chData[data].expGaussFit(ax=ax[1],xLim = self.chData[data].fitwindow)
            else:
                self.chData[data].EnergyHist1DPlot(ax = ax[0,i], ChBins =ChBins, BinRange = Binrange[i], norm = False, log = False,detectors = detectors)#Plots 
                m = self.chData[data].expGaussFit(ax=ax[1,i],xLim = data.fitwindow)

            


            f.write(f'{detectors[f'Channel {self.chData[data].ch}'][0]},{m.values[0]},{m.values[1]},{m.values[2]},{m.values[3]},{m.values[4]},{m.values[5]},{m.values[6]},{m.fval/m.ndof} \n')
        # except:
            #     pass
            
        
        plt.savefig(saveFilePath/f'{fileName}_fit_results.png')
        # plt.show()
        
        f.close()
                
def readInData(file):
    with open(file) as f: 
        for i in range(3): #skips the first 3 lines of the file. 
            next(f)
            
        for i, line in enumerate(f):
            if i == 0: #Look at the third line to get the total number of channels in this coincidence. 
                numChannels = int(line.split("\n")[0].split(": ")[1])
            else:
                data = line.split("\n")[0].split(";")
                ch, E, t,flag = [],[],[],[]
                for j in range(numChannels):
                    ch.append(data[6*j+1])
                    t.append(float(data[6*j+2])/1e3)
                    E.append(float(data[6*j+3]))
                    flag.append(data[6*j+5])
                if i == 1: #Use the first line to initialize the class objects using the channels. 
                    cData = coincData(ch,t,E,flag)
                else:
                    cData.AddEvent(ch,t,E,flag)
                    
    return cData

def BinHistograms(data, Range, numBins):
    #########################################
    #   Small function to bin histograms    #
    #########################################
    
    binRange = np.linspace(Range[0],Range[1],numBins) #Sets the range to bin over. 
    hist,binedges = np.histogram(data,binRange) #Bins the histogram. 
    
    return hist, binedges

def ReadInChannelNames(settings):
    ##############################################################################################
    #   Function that reads in the settings.xml file produced by CoMPASS to grab the channel     #
    #   labels and make a dictionary for them. Note that I did use chatGPT to make this function #
    #   As much as it shames me, I had no idea how to read in this file so I cheated a bit.      #
    #   I have gone through the code, and the comments that ChatGPT left and it all makes sense. #
    #   I have also tested it to ensure that it is working properly.                             #
    ##############################################################################################    

    # Parse the XML file into a tree structure
    tree = ET.parse(settings)
    root = tree.getroot()

    # Initialize an empty dictionary to store channel index and label pairs
    channels = {}

    # Iterate over every <channel> element in the file
    for ch in root.iter("channel"):
        # Find the <index> tag within this <channel> element
        index_elem = ch.find("index")
        
        # Find the <values> tag, which contains parameter entries like labels
        values_elem = ch.find("values")
        
        # Continue only if both index and values are found
        if index_elem is not None and values_elem is not None:
            label = None  # Placeholder for the channel label text

            # Search through each <entry> element under <values>
            for entry in values_elem.findall("entry"):
                # Each entry has a <key> and <value>
                key_elem = entry.find("key")
                value_elem = entry.find("value")
                
                # Check if this entry is the label for the channel
                if key_elem is not None and key_elem.text == "SW_PARAMETER_CH_LABEL":
                    # If found, extract the label text (safely)
                    label = value_elem.text.strip() if value_elem is not None else None
                    break  # No need to check further entries for this channel

            # If a label was found, store it in the dictionary with its index
            if label is not None:
                channels[f'Channel {int(index_elem.text)}'] = [label]

    # Print the resulting dictionary: {channel_index: label, ...}
    return channels


def exp_gauss(x,n_gauss,n_exp,tau,sigma,mu_gauss,mu_exp, xLim1,xLim2):
    return n_gauss+n_exp,(n_gauss*Stats.norm.pdf(x,loc = mu_gauss, scale = sigma) + n_exp * truncexpon.pdf(x, xLim1, xLim2, loc = mu_exp, scale = tau))

def exp_gauss_CDF(x,n_gauss,n_exp,tau,sigma,mu_gauss,mu_exp, xLim1,xLim2):
    return (n_gauss*Stats.norm.cdf(x,loc = mu_gauss, scale = sigma) + n_exp * truncexpon.cdf(x, xLim1, xLim2, loc = mu_exp, scale = tau))

def gauss(x,n_gauss,sigma,mu_gauss):
    return n_gauss*Stats.norm.pdf(x, loc = mu_gauss, scale = sigma)

def expFunc(x, n_exp, mu_exp, tau, xLim1,xLim2):
    return n_exp, n_exp * truncexpon.pdf(x, xLim1, xLim2, loc = mu_exp, scale = tau)

filepath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_06_24/2026_06_24_Daily_LSC_calibration_bck_no_coinc_new_settings_lower_thresh/RAW/coinc_sorted_0ns')
settingsFilePath = filepath.parent.parent / 'settings.xml'
detectors = ReadInChannelNames(settingsFilePath)

# channels = [8,10,12,14]
# fitRegions = [[[2200,3400]],
#               [[1800,3200]],
#               [[1800,3400]],
#               [[1950,3000]]]

# BinRange = [[[0,4000]],
#             [[0,4000]],
#             [[0,4000]],
#             [[0,4000]]]

channels = [2,3,4,5]
fitRegions = [[[3000,4000]],
              [[2500,4000]],
              [[3000,4000]],
              [[2500,4000]]]

BinRange = [[[0,4000]],
            [[0,4000]],
            [[0,4000]],
            [[0,4000]]]
integralBins = 100
scale = False
norm = False

for i,ch in enumerate(channels):
    print(f'Reading in channel: {ch} Data')
    
    if isinstance(ch,(list,tuple)):
        chSuffix = '_'.join(str(c) for c in ch)
    else:
        chSuffix = str(ch)
    
    baseName = filepath.parent.parent.name
    print(baseName)
    fileName = f'SDataR_{baseName}_coinc_{chSuffix}'
    file = filepath / f'{fileName}.txt'

    saveFilePath = filepath / fileName / 'figures'
    saveFilePath.mkdir(parents=True, exist_ok=True)

    
    cdata = readInData(file)
    cdata.Stabilityplots(ChBins = [integralBins,integralBins,integralBins],saveFilePath = saveFilePath, fileName = fileName, scale = scale,detectors=detectors)
    cdata.fitData(ChBins =integralBins , Binrange = BinRange[i], saveFilePath =saveFilePath, fileName=fileName,fitWindow=fitRegions[i], detectors=detectors)