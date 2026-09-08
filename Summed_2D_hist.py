import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import time 
import datetime
from matplotlib import cm
from iminuit import Minuit
from iminuit.cost import ExtendedUnbinnedNLL, ExtendedBinnedNLL
from numba_stats import truncnorm, truncexpon
from scipy.stats import skewnorm
import scipy.stats as Stats


startTime = time.time()

class ChData:
    def __init__(self,ch,E):
        self.ch = int(ch)
        self.E = [E]
        # self.t = [t]
        
    def AddEvent(self,E):
        self.E.append(E)
        # self.t.append(t)
    
        
        
class coincData:
    def __init__(self,channel,Energy):
        self.chList = channel
        self.chData = {}
        for i,ch in enumerate(channel):
            self.chData.update({f'{ch}' : ChData(int(ch),Energy[i])})
        
        # print(self.chData['4'].ch)
        
    def AddEvent(self,channel,Energy):
        
        for i, ch in enumerate(channel):
            self.chData[f'{ch}'].AddEvent(Energy[i])
    

class Summed2DHist:
    def __init__(self,ch):
        self.dates = []
        self.Channels = ch #List of the channels that are being summed (LSC is always along the x-axis. NaI are always along the y-axis)
        self.xInts = np.array([]) #Integrals along the x axis
        self.yInts = np.array([]) #Integrals along the y axis
        
        self.SeparationPoints = []
        self.sumChannels = np.array([])
        
    def addHistData(self,xData,yData, date):

        self.xInts = np.concatenate((self.xInts,xData)) #Adds the latest data to the x and y arrays
        self.yInts = np.concatenate((self.yInts,yData))
        
        
        xData = np.asarray(xData,dtype='float64')
        yData = np.asarray(yData,dtype='float64')
        self.sumChannels = np.concatenate((self.sumChannels,(xData + yData)))
        if date in self.dates:
            self.SeparationPoints[-1] + len(self.xInts)
        else:
            self.dates.append(date)
            self.SeparationPoints.append(len(self.xInts)) #Adds the separation point to note where each day stops. 
        
    
    def plot2DHist(self,savefilepath,scaleData):
        
        bins = 100
        
        if scaleData:
            binRangeX = [0,750]
            binRangeY = [0,750]
        else:
            binRangeX = [0,2250]
            binRangeY = [0,1750]
        fig,ax = plt.subplots(6,6,figsize = (40,40))
        plt.tight_layout()
        plt.subplots_adjust(left = 0.05, bottom = 0.05,hspace = 0.2,wspace = 0.2)
        axsFlat = ax.flatten()
        
        for i,axs in enumerate(axsFlat):
            axs.hist2d(self.xInts[0:self.SeparationPoints[i]],self.yInts[0:self.SeparationPoints[i]], bins = [bins,bins], range = [binRangeX,binRangeY], cmin = 1,rasterized=True,edgecolor='face')
            
            if scaleData:
                axs.set_xlabel('LSC Int (keV)')
                axs.set_ylabel('NaI Int (keV)')
            else:
                axs.set_xlabel('LSC Int (ADC)')
                axs.set_ylabel('NaI Int (ADC)')
            axs.set_title(f'Day {i+1}')
            
        savepath = savefilepath / 'Summed_Histograms' / f'Channels_{self.Channels[0]}_{self.Channels[1]}'
        savepath.mkdir(parents=True,exist_ok=True)
        
        if scaleData:
            plt.savefig(savepath / f'Summed_hist_ch_{self.Channels[0]}_{self.Channels[1]}_scaled_int.png')  
            plt.savefig(savepath / f'Summed_hist_ch_{self.Channels[0]}_{self.Channels[1]}_scaled_int.pdf')  
        else:
            plt.savefig(savepath / f'Summed_hist_ch_{self.Channels[0]}_{self.Channels[1]}.png')  
            plt.savefig(savepath / f'Summed_hist_ch_{self.Channels[0]}_{self.Channels[1]}.pdf')  
        plt.close()
        
    def plotSummed1DHist(self,savefilepath,scaleDat,fitData):
        
        bins = 100
        if scaleData:
            binRange = [0,1400]
        else:
            binRange = [0,4000]
        
        fig,ax = plt.subplots(6,6,figsize = (40,40))
        plt.tight_layout()
        plt.subplots_adjust(left = 0.05, bottom = 0.05,hspace = 0.2,wspace = 0.2)
        axsFlat = ax.flatten()
        

        for i,axs in enumerate(axsFlat):
            counts, mplbins,patches = axs.hist(self.sumChannels[0:self.SeparationPoints[i]], bins = bins, range = binRange,histtype = 'step')
            if fitData:
                self.FitData(axs,binRange,counts,bins)
            if scaleData:
                axs.set_xlabel('LSC Int + NaI Int (keV)')
            else:
                axs.set_xlabel('LSC Int + NaI Int (ADC)')
                
                
            axs.set_ylabel('Counts/bin')
            axs.set_title(f'Day {i+1}')
            
        savepath = savefilepath / 'Summed_Histograms' / f'Channels_{self.Channels[0]}_{self.Channels[1]}'
        savepath.mkdir(parents=True,exist_ok=True)
        if scaleData:
            plt.savefig(savepath / f'Summed_channels_ch_{self.Channels[0]}_{self.Channels[1]}_scaled_int.png')  
            plt.savefig(savepath / f'Summed_channels_ch_{self.Channels[0]}_{self.Channels[1]}_scled_int.pdf')  
        elif fitData:
            plt.savefig(savepath / f'Summed_channels_ch_{self.Channels[0]}_{self.Channels[1]}_fits.png')  
            plt.savefig(savepath / f'Summed_channels_ch_{self.Channels[0]}_{self.Channels[1]}.pdf')  
        else:
            plt.savefig(savepath / f'Summed_channels_ch_{self.Channels[0]}_{self.Channels[1]}.png')  
            plt.savefig(savepath / f'Summed_channels_ch_{self.Channels[0]}_{self.Channels[1]}.pdf')  
        plt.close()
        
    def FitData(self,ax,xLim,counts,nBins):
        
        bins,step = np.linspace(xLim[0],xLim[1],nBins+1,retstep=True)
        
        C = ExtendedBinnedNLL(counts,bins,SkewPlusGauss_CDF)
        
        #1. nSkew
        #2. muSkew
        #3. sigmaSkew
        #4. skew
        #5. nGauss
        #6. muGauss
        #7. sigmaGauss
        self.nSkew = 0.5*sum(counts)
        self.muSkew = 600
        self.sigmaSkew = 0.1*self.muSkew
        self.skew = 4
        self.nGauss = 0.5*sum(counts)
        self.muGauss = 1500
        self.sigmaGauss = 0.1*self.muGauss
        
        initFit = [self.nSkew,self.muSkew,self.sigmaSkew,self.skew,self.nGauss,self.muGauss,self.sigmaGauss]
        
        m = Minuit(C,nSkew = initFit[0],muSkew = initFit[1],sigmaSkew = initFit[2],skew = initFit[3], nGauss = initFit[4], muGauss = initFit[5], sigmaGauss= initFit[6])
        m.limits['nSkew', 'muSkew', 'sigmaSkew', 'nGauss', 'muGauss', 'sigmaGauss'] = (0,None)
        m.limits['skew'] = (10,25)
        
        m.migrad()
        m.hesse()
        
        print(r"$\chi$^2 = " + f"{round(m.fval,2)}")
        print(r"dof = " + f"{round(m.ndof,2)}")
        print(r'$chi$^2 / ndof = ' + f"{m.fval/m.ndof}")
        print(m)
        
        xRange = np.linspace(xLim[0], xLim[1], 1000)
        
        skewPlusGauss = step * SkewPlusGauss_PDF(xRange,m.values[0],m.values[1],m.values[2],m.values[3],m.values[4],m.values[5],m.values[6])
        skewFunc = step * SkewPDF(xRange,m.values[0],m.values[1],m.values[2],m.values[3])
        GaussFunc = step * gauss(xRange,m.values[4],m.values[5],m.values[6]) 
        
        ax.plot(xRange,skewPlusGauss,label = 'fit', linewidth = 4)
        ax.plot(xRange,skewFunc, linestyle = 'dashed', color = 'black', alpha = 0.5)
        ax.plot(xRange,GaussFunc, linestyle = 'dashed', color = 'black', alpha = 0.5)
        # ax.plot([],[],' ', label = r"$\chi$^2 = " + f"{round(m.fval,2)}")
        # ax.plot([],[],' ', label = r"dof = " + f"{round(m.ndof,2)}")
        ax.plot([],[],' ', label = r"$\chi$^2/dof = " + f"{round(m.fval/m.ndof,2)}")    
        
        ax.legend(loc = 'best')
        
        
        return m
        
    def plotOverlayedHist(self,savefilepath,scaleData):
        
        bins = 100
        if scaleData:
            binRange = [0,1400]
        else:
            binRange = [0,4000]
        
        fig,ax = plt.subplots(1,1,figsize = (10,10))
        
        plt.tight_layout()
        plt.subplots_adjust(left = 0.1, bottom = 0.08)
        cmap = cm.viridis
        colours = [cmap(i/(len(self.SeparationPoints))) for i in range(len(self.SeparationPoints))]
        
        for i in range(len(self.SeparationPoints)):
            counts,binsedges = BinHistograms(self.sumChannels[0:self.SeparationPoints[i]],binRange,bins)
            ax.hist(binsedges[:-1],bins = binsedges,weights = counts/max(counts), histtype = 'step', color = colours[i])
            # ax.hist(self.sumChannels[0:self.SeparationPoints[i]], bins = bins, range = binRange, histtype='step', color=colours[i])
        if scaleData:
            ax.set_xlabel('LSC Int + NaI Int (keV)')
        else:
            ax.set_xlabel('LSC Int + NaI Int (ADC)')
        ax.set_ylabel('Counts/bin')
        
        sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=1, vmax=len(self.SeparationPoints)))
        cbar = fig.colorbar(sm, ax=ax)
        cbar.set_label("Day", labelpad=25)
        
        savepath = savefilepath / 'Summed_Histograms' / f'Channels_{self.Channels[0]}_{self.Channels[1]}'
        savepath.mkdir(parents=True,exist_ok=True)
        
        if scaleData:
            plt.savefig(savepath / f'Overlayed_days_ch_{self.Channels[0]}_{self.Channels[1]}_scaled_int.png')
            plt.savefig(savepath / f'Overlayed_days_ch_{self.Channels[0]}_{self.Channels[1]}_scaled_int.pdf')
        else:
            plt.savefig(savepath / f'Overlayed_days_ch_{self.Channels[0]}_{self.Channels[1]}.png')
            plt.savefig(savepath / f'Overlayed_days_ch_{self.Channels[0]}_{self.Channels[1]}.pdf')
        

class ScaleFactors:
    def __init__(self,ch):
        self.ch = ch
        self.scale = {}
        
    def AddData(self,date,scale):
        self.scale[date] = scale
     

def readInData(file, scaleFactor = None):
    with open(file) as f: 
        for i in range(3): #skips the first 3 lines of the file. 
            next(f)
            
        date = str(file.parent.parent.parent.parent).split('/')[-1]
        date = date.replace('_','-')
        for i, line in enumerate(f):
            if i == 0: #Look at the third line to get the total number of channels in this coincidence. 
                numChannels = int(line.split("\n")[0].split(": ")[1])
            else:
                
                data = line.split("\n")[0].split(";")
                ch, E, t = [],[],[]
                for j in range(numChannels):
                    ch.append(int(data[6*j+1]))
                    # t.append(float(data[6*j+2])/1e3)
                    if scaleFactor is not None:
                        for s in scaleFactor:
                            if s.ch == int(data[6*j+1]):
                                scaleFac = s.scale[date]
                                break
                            else:
                                pass
                                
                        E.append(float(data[6*j+3])/scaleFac)
                        
                    else:
                        E.append(float(data[6*j+3]))
                if i == 1: #Use the first line to initialize the class objects using the channels. 
                    cData = coincData(ch,E)
                else:
                    cData.AddEvent(ch,E)
                    
    return cData



def BinHistograms(data, Range, numBins):
    #########################################
    #   Small function to bin histograms    #
    #########################################
    
    binRange = np.linspace(Range[0],Range[1],numBins) #Sets the range to bin over. 
    hist,binedges = np.histogram(data,binRange) #Bins the histogram. 
    
    return hist, binedges

def readInScaleFactors(filepath, channels):
    
    scale = []
    for ch in channels:
        scale.append(ScaleFactors(ch))
    with open(filepath) as f:
        next(f) #Skips the header
        
        for line in f:
            data = line.split(',')
            for s in scale:
                if s.ch == int(data[1]):
                    s.AddData(data[0].split(' 00:00:00')[0], float(data[3]))

    return scale


def SkewPlusGauss_CDF(x,nSkew,muSkew,sigmaSkew,skew,nGauss,muGauss,sigmaGauss):
    skew = nSkew*skewnorm.cdf(x,a = skew, loc = muSkew, scale = sigmaSkew)
    Gauss = nGauss*Stats.norm.cdf(x, loc = muGauss, scale = sigmaGauss)
    
    return skew + Gauss

def SkewPlusGauss_PDF(x,nSkew,muSkew,sigmaSkew,skew,nGauss,muGauss,sigmaGauss):
    skew = nSkew*skewnorm.pdf(x,a = skew, loc = muSkew, scale = sigmaSkew)
    Gauss = nGauss*Stats.norm.pdf(x, loc = muGauss, scale = sigmaGauss)
    
    return skew + Gauss

def SkewPDF(x,nSkew,muSkew,sigmaSkew,skew):
    return nSkew*skewnorm.pdf(x,a = skew, loc = muSkew, scale = sigmaSkew)  

def gauss(x,nGauss,muGauss,sigmaGauss):
    return nGauss*Stats.norm.pdf(x, loc = muGauss, scale = sigmaGauss)
            # scaleFactors[f'{data[0]}'] = []
# def readInData(file):
#     with open(file) as f:
#         lines = f.readlines()

#     numChannels = int(lines[3].split(": ")[1])

#     # Parse all data rows at once
#     data_lines = [line.strip().split(";") for line in lines[4:] if line.strip()]
#     raw = np.array(data_lines)

#     # Extract all channels, times, energies as vectors
#     ch_cols = [raw[:, 6*j + 1] for j in range(numChannels)]
#     t_cols  = [raw[:, 6*j + 2].astype(float) / 1e3 for j in range(numChannels)]
#     E_cols  = [raw[:, 6*j + 3].astype(float) for j in range(numChannels)]

#     channels = [int(col[0]) for col in ch_cols]

#     # Initialize with first event (preserves your existing constructor signature)
#     cData = coincData(channels, [t[0] for t in t_cols], [E[0] for E in E_cols])

#     # Bulk-load remaining events into ChData, bypassing the slow append loop
#     for ch, E_vec, t_vec in zip(channels, E_cols, t_cols):
#         cData.chData[ch].E = list(E_vec)
#         cData.chData[ch].t = list(t_vec)

#     return cData

    
    
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


rootfilePath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/') #filepath to the coinc sorted directory. 

#A list of data sets that are excluded from the data. This can be due to bad data or incorrect settings.
excludedDataFiles = [Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/2026_04_28/2026_04_28_Daily_LSC_calibration_Cs137_coinc'), #Data set was taken before the settings were finalized
                     ]

averageScaleFactorFP = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/Results/Annulus_stability_data_average.txt')

LSCChannels = [4,5]
NaIChannels = [8,10,12,14] #Protects against the possibility of having a werid channel coincidence layout with incorrect channel numbers. 

# LSCChannels = [0,1]
# NaIChannels = [2,3,4,5]

scaleData = True
fitData = False
if scaleData:
    scaleFac = readInScaleFactors(averageScaleFactorFP, channels = LSCChannels + NaIChannels)

pattern = re.compile(rf'_coinc_{LSCChannels[0]}_{LSCChannels[1]}_({NaIChannels[0]}|{NaIChannels[1]}|{NaIChannels[2]}|{NaIChannels[3]}).txt$')

coincFiles = sorted([
    f for f in rootfilePath.glob(f'**/*_coinc_{LSCChannels[0]}_{LSCChannels[1]}_*.txt')
    if pattern.search(f.name)
])
histogramData = []
for i in LSCChannels:
    for j in NaIChannels:
        histogramData.append(Summed2DHist([i,j]))
        print(histogramData[-1].Channels)


readInStartTime = time.time()
print('Reading in histogram data:')

for i,filePath in enumerate(coincFiles):
    if filePath.parent.parent.parent in excludedDataFiles:
        print(f'Skipped file: {filePath}')
        pass
    else:
        date = filePath.parent.parent.parent.parent.stem
        settingsFilePath = filePath.parent.parent.parent / 'settings.xml'
        # saveFilePath = filePath.parent.parent / 'Daily_calibration_fits' / 'figures'
        # outputFileName = filePath.stem
        # outputFile = saveFilePath / f"{outputFileName}_Profile_hist_fit_results.txt"
        
        # if outputFile.exists():
        #     outputFile.unlink()

        # saveFilePath.mkdir(parents=True, exist_ok=True)
        savefilepath = Path('/home/nick/PhD/KDK+/Daily_LSC_Calibration_testing/stability_figures/')
        if i == 0:
            Detectors = ReadInChannelNames(settingsFilePath)
        
        if scaleData:
            cData = readInData(filePath,scaleFac)
        else:
            cData = readInData(filePath)
            
        chPairs = [[cData.chList[0],cData.chList[2]], [cData.chList[1],cData.chList[2]]] #Makes a 2D array with the two channel pairs for the summed histogram data. 
        print(cData.chList)
        # print(len(cData.chData['4'].E))


        for j,hist in enumerate(histogramData):
            if hist.Channels == chPairs[0]:
                hist.addHistData(cData.chData[str(chPairs[0][0])].E,cData.chData[str(chPairs[0][1])].E,date)
            elif hist.Channels == chPairs[1]:
                hist.addHistData(cData.chData[str(chPairs[1][0])].E,cData.chData[str(chPairs[1][1])].E,date)


totalReadTime = time.time()
print(f'Total time to read in data: \t {totalReadTime - readInStartTime} s')

for hist in histogramData:
    print(f'Plotting Histogram from channels: {hist.Channels}')
    hist.plot2DHist(savefilepath,scaleData)
    print(f"Plotting Summed Channel data from channels: {hist.Channels}")
    hist.plotSummed1DHist(savefilepath,scaleData,fitData)
    print(f"Plotting overlay data for channels: {hist.Channels}")
    hist.plotOverlayedHist(savefilepath,scaleData)



    
    
    
    