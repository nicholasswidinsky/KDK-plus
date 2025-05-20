import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path

matplotlib.rcParams["font.size"] = 15
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["mathtext.default"] = 'regular'
matplotlib.rcParams['lines.markersize'] = 3



def BinHistograms(data, Range, numBins):
    binRange = np.linspace(Range[0],Range[1],numBins)
    hist,binedges = np.histogram(data,binRange)
    
    return hist, binedges


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
                    data[c][0].append(float(lines[2])/(1e12))
                    data[c][1].append(float(lines[3]))
                    data[c][2].append(lines[5])
    return data, Channels


def TimeDifferenceChannel(data,channels):
    ######################################################################################################
    # Data: 2D array that contains the data of all the channels pre sorted using the RedInFile function. #
    # Channels: array of the two channels that you are taking the time difference of.                    #
    ######################################################################################################
    
    timeDiff = []
    # print(len(data[channels[0]][0]))
    # print(len(data[channels[1]][0]))
    
    for i in range(len(data[0][0])):
        timeDiff.append(data[channels[0]][0][i] - data[channels[1]][0][i])
        
    timeDiff = np.asfarray(timeDiff, float)
    return timeDiff



def energyHist1D(data, Bckdata, ChBins, BinRange, norm, log,saveFilePath,fileName,Channel, bck):
    ######################################################################################################
    # Data: 2D array that contains the data of all the channels pre sorted using the RedInFile function. #
    # bckData: 2D array that contains the data of the background run.                                    #
    # ChBins: Array of the number of bins you want for each channel                                      #
    # BinRange: Range that the histogram will be binned over.                                            #
    ######################################################################################################
    # bins = 100
    # binRange = [0,1000]

    # binsLYSO = 100
    # binRangeLYSO = [0,2000]
    
    root = np.sqrt(len(data))
    if root == int(root):
        nrows = root
        ncols = root
    
    else: 
        nrows = int(root)
        ncols = len(data)//nrows 
        
        if (ncols * nrows) != len(data):
            ncols +=1
        
        
        
    fig, ax = plt.subplots(nrows,ncols, figsize = (10*ncols,10*nrows))
    fig.tight_layout()
    plt.subplots_adjust(wspace = 0.1,hspace =0.3)
    
    if nrows == 1:
        for i in range(len(data)):
            Int, Bins = BinHistograms(data[i][1],BinRange[i],ChBins[1])

            ax[i].hist(Bins[:-1],bins =Bins,density = norm,weights=Int/(data[i][0][-1] - data[i][0][0]),histtype = 'step',label = f'Ch {Channel[i]}', log = log)
        
            if bck:
                bckInt, bckBins = BinHistograms(Bckdata[i][1],BinRange[i],ChBins[1])
                ax[i].hist(bckBins[:-1],bins =bckBins,density = norm,weights=bckInt/(Bckdata[i][0][-1] - Bckdata[i][0][0]),histtype = 'step',label = f'Background Ch {Channel[i]} ', log = log)

            ax[i].set_title(f'Channel {Channel[i]}')
            ax[i].set_xlabel('Integral (ADC Channel)')
            ax[i].set_ylabel('Counts/bin/ps')
            ax[i].legend(loc = 'best')

        
    else: 
        row = 0
        col = 0
        
        for i in range(len(data)):
            Int, Bins = BinHistograms(data[i][1],BinRange[i],ChBins[1])

            ax[row,col].hist(Bins[:-1],bins =Bins,density = norm,weights=Int/(data[i][0][-1] - data[i][0][0]),histtype = 'step',label = f'Ch {Channel[i]}', log = log)
            
            
            if bck:
                bckInt, bckBins = BinHistograms(Bckdata[i][1],BinRange[i],ChBins[1])
                ax[row,col].hist(bckBins[:-1],bins =bckBins,density = norm,weights=bckInt/(Bckdata[i][0][-1] - Bckdata[i][0][0]),histtype = 'step',label = f'Background Ch {Channel[i]} ', log = log)

            ax[row,col].set_title(f'Channel {Channel[i]}')
            ax[row,col].set_xlabel('Integral (ADC Channel)')
            ax[row,col].set_ylabel('Counts/bin/ps')
            
            row += 1
            if row == nrows:
                row = 0
                col +=1
        
    plt.savefig(f'{saveFilePath}/{fileName}_integral_1D_hist.png',bbox_inches='tight')
    plt.close()
    
def Stability_plots(data, ChBins, BinRange,saveFilePath,fileName,Channel):
    
    
    
    fig, ax = plt.subplots(2,len(Channels), figsize = (10*len(Channels),20))
    plt.tight_layout()
    plt.subplots_adjust(left = 0.06,wspace = 0.15,hspace = 0.1,top = 0.98,bottom = 0.04)
    for i in range(len(Channels)):
        energybinsRange = np.linspace(BinRange[i][0],BinRange[i][1],ChBins[0])
        

        timebinsRange = np.linspace(min(data[i][0]),max(data[i][0]),100)
        
        ax[0,i].hist(data[i][1],bins = energybinsRange)
        ax[0,i].set_title(f'Channel {Channel[i]}')
        ax[0,i].set_xlabel('Integral (ADC Channel)')
        ax[0,i].set_ylabel('Counts/bin')
        
        
        ax[1,i].hist2d(data[i][0],data[i][1], bins = (timebinsRange,energybinsRange))
        ax[1,i].set_title(f'Channel {Channel[i]}')
        ax[1,i].set_xlabel('Trigger Time (s)')
        
        ax[1,i].set_ylabel('Integral (ADC Channel)')
    plt.savefig(f'{saveFilePath}/{fileName}_stability_plots.png',bbox_inches='tight')
    
    
def TimeDiff1D(data, BinRange,saveFilePath,fileName,Channels):
    ######################################################################################################
    # Data: 2D array that contains the data of all the channels pre sorted using the RedInFile function. #
    # bckData: 2D array that contains the data of the background run.                                    #
    # timeBins: Array of the number of bins you want for each channel and time bin                       #
    # BinRange: Range that the histogram will be binned over.  
    ######################################################################################################   
    
    
    
    timeBins = (BinRange[1] - BinRange[0])//2000

    timebinsRange = np.linspace(BinRange[0],BinRange[1],timeBins)/1000


    fig,ax = plt.subplots(len(data),len(data),figsize = (10*len(data),10*len(data)))
    
    for i in range(len(data)):
        for j in range(len(data)):

            ax[i,j].hist(TimeDifferenceChannel(data,[i,j])/1000, bins = timebinsRange, log = True)
            ax[i,j].set_title(f'Time Difference Ch{Channels[i]} - Ch{Channels[j]}')
            ax[i,j].set_xlabel('Time Difference (ns)')
            ax[i,j].set_ylabel('Counts')

    


    # plt.show()
    plt.savefig(f'{saveFilePath}/{fileName}_time_diff_1D_hist.png',bbox_inches='tight')
    plt.close()
    
def energyHist2D(data,ChBins,BinRange,saveFilePath,fileName,Channels,log):
    
    ChBinsArray = []
    for i in range(len(data)):
        ChBinsArray.append(np.linspace(BinRange[i][0],BinRange[i][1],ChBins[i]))
    
    
    
    fig,ax = plt.subplots(len(data),len(data),figsize = (10*len(data),10*len(data)))    
    
    for i in range(len(data)):
        for j in range(len(data)):
            
            if log:
                ax[i,j].hist2d(data[i][1],data[j][1],bins = (ChBinsArray[i],ChBinsArray[j]),cmin = 1,norm=matplotlib.colors.LogNorm())
            else:
                ax[i,j].hist2d(data[i][1],data[j][1],bins = (ChBinsArray[i],ChBinsArray[j]),cmin = 1)
            ax[i,j].set_title(f'Ch{Channels[i]}, Ch{Channels[j]} 2D hist')
            ax[i,j].set_xlabel(f'Ch{Channels[i]} Integral')
            ax[i,j].set_ylabel(f'Ch{Channels[j]} Integral')
    plt.savefig(f'{saveFilePath}/{fileName}_integral_2D_hist.png',bbox_inches='tight')
    plt.close()
    
def TimeDiff2D(data,ChBins,ChBinRange,TimeBinRange,saveFilePath,fileName,Channels,log):
    
    ChBinsArray = []
    for i in range(len(data)):
        ChBinsArray.append(np.linspace(ChBinRange[i][0],ChBinRange[i][1],ChBins[i]))
        
    timeBins = (TimeBinRange[1] - TimeBinRange[0])//2000
    timebinsRange = np.linspace(TimeBinRange[0],TimeBinRange[1],timeBins)/1000
    
    
    fig,ax=plt.subplots(len(data),len(data),figsize = (10*len(data),10*len(data)))
    
    for i in range(len(data)):
        for j in range(len(data)):
            
            if log:
                ax[i,j].hist2d(data[i][1],TimeDifferenceChannel(data,[i,j])/1000, bins =(ChBinsArray[i],timebinsRange), cmin = 1,norm=matplotlib.colors.LogNorm())
            else:
                ax[i,j].hist2d(data[i][1],TimeDifferenceChannel(data,[i,j])/1000, bins =(ChBinsArray[i],timebinsRange), cmin = 1)
            ax[i,j].set_title(f'Ch {Channels[i]} Integral, Ch{Channels[i]} - Ch{Channels[j]} Time Difference')
            ax[i,j].set_xlabel(f'Ch{Channels[i]} integral')
            ax[i,j].set_ylabel(f'Ch{Channels[i]} - Ch{Channels[j]} Time Difference')      
    plt.savefig(f'{saveFilePath}/{fileName}_integral_time_diff_2D_hist.png',bbox_inches='tight')
    plt.close()
    
                
        

def TimeCut(data,timeCut,CoincCH,Channels):
    ######################################################################################################
    # Data: N-dimensional array that contains three more arrays per element. These sub arrays hold the information for each event. The N-dimensions of the total array correspond to individual channels. 
    # 
    # timeCut: array with two elements that hold the start and the end of the timecut. 
    #
    # CoincCH: Array with two elements that corresponds to the two channels that will be compared for the time cut. 
    #
    # Channels: Array of all the channels used in the experiment. This is done purely to coordinate the indeces for data to the CoincCH values. 
    ######################################################################################################
    
    Channels = np.asfarray(Channels, int)
    CoincCh = np.asfarray(CoincCh, int)

    Channelind = np.where((Channels >= sorted(CoincCh)[0]) & (Channels <= sorted(CoincCh)[1]))
    
    timeDiff = TimeDifferenceChannel(data,Channelind)
    
    CoincData = [[],[],[]]
    cutData = [[],[],[]]
    for i in range(len(data[0][0])):
        if timeDiff[i] > timeCut[0] and timeDiff[i] < timeCut[1]:
            CoincData[0].append(data[0][1][i])
            CoincData[1].append(data[1][1][i])
            CoincData[2].append(data[2][1][i])
        else:
            cutData[0].append(data[0][1][i])
            cutData[1].append(data[1][1][i])    
            cutData[2].append(data[2][1][i])    
            
    
    CoincData[0] = np.asfarray(CoincData[0],float)
    CoincData[1] = np.asfarray(CoincData[1],float)
    CoincData[2] = np.asfarray(CoincData[2],float)
    
    cutData[0] = np.asfarray(cutData[0],float)
    cutData[1] = np.asfarray(cutData[1],float)
    cutData[2] = np.asfarray(cutData[2],float)
    
    # CoincData = np.asfarray(CoincData,float)
    # cutData = np.asfarray(cutData,float)
    
    return CoincData, cutData
                

        
    
def TimeCutHist(data,ChBins,ChBinRange,norm,TimeBinRange,timeCut,saveFilePath,fileName):
    
    CoincData, cutData = TimeCut(data,timeCut)
    # print(type(CoincData))
    # print(CoincData)
    
    fig,ax = plt.subplots(1,3, figsize = (30,10))
    
    # Int0, Bins0 = BinHistograms(data[0][1],BinRange[0],ChBins[0])
    # Int2, Bins2 = BinHistograms(data[1][1],BinRange[1],ChBins[1])
    # Int4, Bins4 = BinHistograms(data[2][1],BinRange[2],ChBins[2])

    # bckInt0, bckBins0 = BinHistograms(Bckdata[0][1],BinRange[0],ChBins[0])
    # bckInt2, bckBins2 = BinHistograms(Bckdata[1][1],BinRange[1],ChBins[1])
    # bckInt4, bckBins4 = BinHistograms(Bckdata[2][1],BinRange[2],ChBins[2])
    
    # ax[0].hist(Bins0[:-1],bins =Bins0,density = normalize,weights=Int0/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Ch 0 (LSC)', log = logScale)
    
    for i in range(3):
        Int,Bins = BinHistograms(CoincData[i], ChBinRange[i],ChBins[i])
        ax[i].hist(Bins[:-1],bins = Bins,density = norm,weights = Int/(data[i][0][-1]-data[i][0][0]),histtype = 'step', label = 'Time Cut', log = True)
        
        IntRaw,BinsRaw = BinHistograms(data[i][1], ChBinRange[i],ChBins[i])
        ax[i].hist(BinsRaw[:-1],bins = BinsRaw,density = norm,weights = IntRaw/(data[i][0][-1]-data[i][0][0]),histtype = 'step', label = 'No Time Cut')
        
        IntNoCut,BinsNoCut = BinHistograms(cutData[i], ChBinRange[i],ChBins[i])
        ax[i].hist(BinsNoCut[:-1],bins = BinsNoCut,density = norm,weights = IntNoCut/(data[i][0][-1]-data[i][0][0]),histtype = 'step', label = 'Cut Data')
        
        ax[i].set_xlabel('Integral (ADC Channel)')
        ax[i].set_ylabel('Counts/bin/ps')
        ax[i].legend(loc = 'best')
    
    ax[0].set_title('Left PMT')
    ax[1].set_title('Right PMT')
    ax[2].set_title('LYSO')
    
    
    plt.savefig(f'{saveFilePath}/{fileName}_time_cut_hist.png',bbox_inches='tight')
    

    timeBins = (TimeBinRange[1] - TimeBinRange[0])//2000
    timebinsRange = np.linspace(TimeBinRange[0],TimeBinRange[1],timeBins)/1000
    
    fig,ax = plt.subplots(3,3,figsize = (30,30))
    
    for i in range(3):
        for j in range(3):

            ax[i,j].hist(TimeDifferenceChannel(data,[i,j])/1000, bins = timebinsRange, log = True)
            ax[i,j].set_title(f'Time Difference Ch{2*i} - Ch{2*j}')
            
            ax[i,j].axvline(timeCut[0]/1000, color = 'red')
            ax[i,j].axvline(timeCut[1]/1000, color = 'red')
            
            
            ax[i,j].set_xlabel('Time Difference (ns)')
            ax[i,j].set_ylabel('Counts')
            
    plt.savefig(f'{saveFilePath}/{fileName}_time_hist_with_cuts.png',bbox_inches='tight')
    
    
    
filepath = '/home/nick/PhD/Plastic_scintillators/Cryostat_gamma_V2/Darkbox_tests/CAEN_DAQ/2025_05_16/'

GammaDetector = 'NaI'
file = "SDataR_Small_PSC_Darkbox_NaI_Cs137_coinc_10CG_PSC_Lead_no_brass.CSV"
bckfile = '/home/nick/PhD/Plastic_scintillators/Cryostat_gamma_V2/Darkbox_tests/CAEN_DAQ/2025_05_13/SDataR_Small_PSC_Darkbox_NaI_bck_coinc_2.CSV'
# bckfile = '/home/nick/PhD/Plastic_scintillators/Cryostat_gamma_V2/Darkbox_tests/CAEN_DAQ/2025_05_15/SDataR_Small_PSC_Darkbox_NaI_Cs137_coinc_10CG_PSC_Lead.CSV'
# bckfile = f'{filepath}/{file}'

# Channels = [0,2]

# fileData = ReadInFileNaI(f'{filepath}/{file}', 2, Channels = Channels)
dataFile, Channels = ReadInFileNaIChannels(f'{filepath}/{file}')

bckdataFile, bckChannels = ReadInFileNaIChannels(f'{bckfile}')


backgrounds = True
normalize = False
Log = False


fileName = file.split('.CSV')[0]
bckfileName = bckfile.split('/')[-1].split('.CSV')[0]

saveFilePath = f"{filepath}/{fileName}/figures"
Path(f"{saveFilePath}").mkdir(parents=True, exist_ok=True)


energyHist1D(data = dataFile, Bckdata = bckdataFile, ChBins = [100,100], BinRange = [[0,4000],[0,1500]], norm = normalize, log = Log,saveFilePath = saveFilePath,fileName = fileName, Channel = Channels, bck = backgrounds)

Stability_plots(data = dataFile, ChBins = [100,100], BinRange = [[0,3000],[0,1000]],saveFilePath = saveFilePath,fileName = fileName, Channel = Channels)
try:
    TimeDiff1D(data = dataFile, BinRange = [-500000,500000],saveFilePath = saveFilePath,fileName = fileName,Channels = Channels)
    energyHist2D(dataFile, ChBins = [100,100], BinRange = [[0,5000],[0,1500]],saveFilePath = saveFilePath,fileName = fileName,Channels = Channels,log = Log)
    TimeDiff2D(data = dataFile,ChBins = [100,100],ChBinRange= [[0,1000],[0,1500]],TimeBinRange= [-400000,400000],saveFilePath = saveFilePath,fileName = fileName, Channels = Channels, log = Log)
except:
    pass

# energyHist2D(bckdataFile, ChBins = [100,100,100,100], BinRange = [[0,4000],[0,4000],[0,4000],[0,4000]],saveFilePath = saveFilePath,fileName = bckfileName, Channels = [0,2,6,8])
# TimDiff2D(data = bckdataFile,ChBins = [100,100,100,100],ChBinRange= [[0,4000],[0,4000],[0,4000],[0,4000]],TimeBinRange= [-400000,400000],saveFilePath = saveFilePath,fileName = bckfileName)


# TimeCutHist(data=dataFile,ChBins = [100,100,100],ChBinRange = [[0,4000],[0,4000],[0,4000]],norm = False,TimeBinRange= [200000,350000],timeCut = [230000,250000],saveFilePath = saveFilePath,fileName= fileName)
