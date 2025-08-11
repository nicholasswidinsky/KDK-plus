import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path

matplotlib.rcParams["font.size"] = 15
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["mathtext.default"] = 'regular'
matplotlib.rcParams['lines.markersize'] = 3

def ReadInFile(file, numChannel):
    
    data = []
    for c in range(numChannel):
        data.append([ [] for _ in range(3)])
    channels = [2*n for n in range(numChannel)]
    
    
    with open(file) as f:
        next(f)
        for line in f:
            lines = line.split('\n')[0].split(';')
            
            for c in range(numChannel):
                if int(lines[1]) == channels[c]:
                    data[c][0].append(float(lines[2]))
                    data[c][1].append(float(lines[3]))
                    data[c][2].append(lines[5])
    return data

def ReadInFileNaI(file, numChannel):
    
    data = []
    for c in range(numChannel):
        data.append([ [] for _ in range(3)])
    # channels = [2*n for n in range(numChannel)]
    channels = [0,2,6]
    
    
    with open(file) as f:
        next(f)
        for line in f:
            lines = line.split('\n')[0].split(';')
            
            for c in range(numChannel):
                if int(lines[1]) == channels[c]:
                    data[c][0].append(float(lines[2]))
                    data[c][1].append(float(lines[3]))
                    data[c][2].append(lines[5])
    return data
    
    
def BinHistograms(data, Range, numBins):
    binRange = np.linspace(Range[0],Range[1],numBins)
    hist,binedges = np.histogram(data,binRange)
    
    return hist, binedges
    
def TimeDifference(data):
    timeDiff0_2 =[]
    
    for i in range(len(data[0][0])):
        timeDiff0_2.append(data[0][0][i]-data[1][0][i])


        
    timeDiff0_2 = np.asfarray(timeDiff0_2,float)


        
    return timeDiff0_2

def BinHistograms(data, Range, numBins):
    binRange = np.linspace(Range[0],Range[1],numBins)
    hist,binedges = np.histogram(data,binRange)
    
    return hist, binedges


def TimeDifferenceChannel(data,channels):
    ######################################################################################################
    # Data: 2D array that contains the data of all the channels pre sorted using the RedInFile function. #
    # Channels: array of the two channels that you are taking the time difference of.                    #
    ######################################################################################################
    
    timeDiff = []
    for i in range(len(data[0][0])):
        timeDiff.append(data[channels[0]][0][i] - data[channels[1]][0][i])
        
    timeDiff = np.asfarray(timeDiff, float)
    return timeDiff

# def energyHist1D(data, Bckdata, ChBins, BinRange, norm, log,saveFilePath,fileName, bck):
#     ######################################################################################################
#     # Data: 2D array that contains the data of all the channels pre sorted using the RedInFile function. #
#     # bckData: 2D array that contains the data of the background run.                                    #
#     # ChBins: Array of the number of bins you want for each channel                                      #
#     # BinRange: Range that the histogram will be binned over.                                            #
#     ######################################################################################################
#     # bins = 100
#     # binRange = [0,1000]

#     # binsLYSO = 100
#     # binRangeLYSO = [0,2000]


#     Int0, Bins0 = BinHistograms(data[0][1],BinRange[0],ChBins[0])
#     Int2, Bins2 = BinHistograms(data[1][1],BinRange[1],ChBins[1])
#     Int4, Bins4 = BinHistograms(data[2][1],BinRange[2],ChBins[2])

#     if bck:
#         bckInt0, bckBins0 = BinHistograms(Bckdata[0][1],BinRange[0],ChBins[0])
#         bckInt2, bckBins2 = BinHistograms(Bckdata[1][1],BinRange[1],ChBins[1])
#         bckInt4, bckBins4 = BinHistograms(Bckdata[2][1],BinRange[2],ChBins[2])
    
#     normalize = norm
#     logScale = log
    
#     fig, ax = plt.subplots(1,3, figsize = (30,10))

#     ax[0].hist(Bins0[:-1],bins =Bins0,density = normalize,weights=Int0/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Ch 0 (LSC)', log = logScale)
#     if bck:
#         ax[0].hist(bckBins0[:-1],bins =bckBins0,density = normalize,weights=bckInt0/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = 'Background Ch 0 (LSC)', log = logScale)

#     ax[1].hist(Bins2[:-1],bins =Bins2,density = normalize,weights=Int2/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Ch 2 (LSC)', log = logScale)
#     if bck:
#         ax[1].hist(bckBins2[:-1],bins =bckBins2,density = normalize,weights=bckInt2/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = 'Background Ch 0 (LSC)', log = logScale)

#     ax[2].hist(Bins4[:-1],bins =Bins4,density = normalize,weights=Int4/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Ch 4 (LYSO)', log = logScale)
#     if bck:
#         ax[2].hist(bckBins4[:-1],bins =bckBins4,density = normalize,weights=bckInt4/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = 'Background Ch 0 (LSC)', log = logScale)

#     ax[0].set_title('Left PMT')
#     ax[1].set_title('Right PMT')
#     ax[2].set_title('LYSO')

#     ax[0].set_xlabel('Integral (ADC Channel)')
#     ax[1].set_xlabel('Integral (ADC Channel)')
#     ax[2].set_xlabel('Integral (ADC Channel)')

#     ax[0].set_ylabel('Counts/bin/ps')
#     ax[1].set_ylabel('Counts/bin/ps')
#     ax[2].set_ylabel('Counts/bin/ps')
    
#     ax[0].legend(loc = 'best')
#     ax[1].legend(loc = 'best')
#     ax[2].legend(loc = 'best')
    
    
#     # plt.show()
#     plt.savefig(f'{saveFilePath}/{fileName}_integral_1D_hist.png',bbox_inches='tight')
#     plt.close()
    
def TimeDiff1D(data, BinRange,saveFilePath,fileName):
    ######################################################################################################
    # Data: 2D array that contains the data of all the channels pre sorted using the RedInFile function. #
    # bckData: 2D array that contains the data of the background run.                                    #
    # timeBins: Array of the number of bins you want for each channel and time bin                       #
    # BinRange: Range that the histogram will be binned over.  
    ######################################################################################################   
    
    
    
    timeBins = (BinRange[1] - BinRange[0])//2000

    timebinsRange = np.linspace(BinRange[0],BinRange[1],timeBins)/1000


    fig,ax = plt.subplots(3,3,figsize = (30,30))
    
    for i in range(3):
        for j in range(3):

            ax[i,j].hist(TimeDifferenceChannel(data,[i,j])/1000, bins = timebinsRange, log = True)
            ax[i,j].set_title(f'Time Difference Ch{2*i} - Ch{2*j}')

            
    ax[0,0].set_xlabel('Time Difference (ns)')
    ax[0,1].set_xlabel('Time Difference (ns)')
    ax[0,2].set_xlabel('Time Difference (ns)')
    ax[1,0].set_xlabel('Time Difference (ns)')
    ax[1,1].set_xlabel('Time Difference (ns)')
    ax[1,2].set_xlabel('Time Difference (ns)')
    ax[2,0].set_xlabel('Time Difference (ns)')
    ax[2,1].set_xlabel('Time Difference (ns)')
    ax[2,2].set_xlabel('Time Difference (ns)')
    
    ax[0,0].set_ylabel('Counts')
    ax[0,1].set_ylabel('Counts')
    ax[0,2].set_ylabel('Counts')
    ax[1,0].set_ylabel('Counts')
    ax[1,1].set_ylabel('Counts')
    ax[1,2].set_ylabel('Counts')
    ax[2,0].set_ylabel('Counts')
    ax[2,1].set_ylabel('Counts')
    ax[2,2].set_ylabel('Counts')
    


    # plt.show()
    plt.savefig(f'{saveFilePath}/{fileName}_time_diff_1D_hist.png',bbox_inches='tight')
    plt.close()
    
def energyHist2D(data,ChBins,BinRange,saveFilePath,fileName):
    
    ChBinsArray = []
    for i in range(3):
        ChBinsArray.append(np.linspace(BinRange[i][0],BinRange[i][1],ChBins[i]))
    
    
    
    fig,ax = plt.subplots(3,3,figsize = (30,30))    
    
    for i in range(3):
        for j in range(3):

            ax[i,j].hist2d(data[i][1],data[j][1],bins = (ChBinsArray[i],ChBinsArray[j]),cmin = 1)#,norm=matplotlib.colors.LogNorm())
            if i == 0 and j == 1:
                ax[i,j].plot(data[i][1],data[i][1],linestyle = "--", c = 'black')
            if i == 1 and j == 0:
                ax[i,j].plot(data[i][1],data[i][1],linestyle = "--", c = 'black')
                
            ax[i,j].set_title(f'Ch{2*i}, Ch{2*j} 2D hist')
            ax[i,j].set_xlabel(f'Ch{2*i} Integral')
            ax[i,j].set_ylabel(f'Ch{2*j} Integral')

    # plt.show()
    plt.savefig(f'{saveFilePath}/{fileName}_integral_2D_hist.png',bbox_inches='tight')
    plt.close()
    
def TimDiff2D(data,ChBins,ChBinRange,TimeBinRange,saveFilePath,fileName):
    
    ChBinsArray = []
    for i in range(3):
        ChBinsArray.append(np.linspace(ChBinRange[i][0],ChBinRange[i][1],ChBins[i]))
        
    timeBins = (TimeBinRange[1] - TimeBinRange[0])//2000
    timebinsRange = np.linspace(TimeBinRange[0],TimeBinRange[1],timeBins)/1000
    
    
    fig,ax=plt.subplots(3,3,figsize = (30,30))
    
    for i in range(3):
        for j in range(3):

            ax[i,j].hist2d(data[i][1],TimeDifferenceChannel(data,[i,j])/1000, bins =(ChBinsArray[i],timebinsRange), cmin = 1,norm=matplotlib.colors.LogNorm())
            ax[i,j].set_title(f'Ch{2*i} Integral, Ch{2*i} - Ch{2*j} Time Difference')
            ax[i,j].set_xlabel(f'Ch{2*i} integral')
            ax[i,j].set_ylabel(f'Ch{2*i} - Ch{2*j} Time Difference')

            
    # plt.show()
    plt.savefig(f'{saveFilePath}/{fileName}_integral_time_diff_2D_hist.png',bbox_inches='tight')
    plt.close()
    
def TimeCut(data):
    timeDiff = TimeDifferenceChannel(data,[4,0])
    
    CoincData = [[],[],[]]
    cutData = [[],[],[]]
    for i in range(len(data)):
        if timeDiff > 50 and timeDiff < 100:
            CoincData[0].append(data[0][i])
            CoincData[1].append(data[1][i])
            CoincData[2].append(data[2][i])
        else:
            cutData[0].append(data[0][i])
            cutData[1].append(data[1][i])    
            cutData[2].append(data[2][i])    
    return CoincData, cutData
                
        
def TimeCut(data,timeCut):
    timeDiff = TimeDifferenceChannel(data,[2,0])
    
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
                
def energyHist1D(data, Bckdata, ChBins, BinRange, norm, log,saveFilePath,fileName,gammaDet,bck):
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


    Int0, Bins0 = BinHistograms(data[0][1],BinRange[0],ChBins[0])
    Int2, Bins2 = BinHistograms(data[1][1],BinRange[1],ChBins[1])
    Int4, Bins4 = BinHistograms(data[2][1],BinRange[2],ChBins[2])

    if bck:
        bckInt0, bckBins0 = BinHistograms(Bckdata[0][1],BinRange[0],ChBins[0])
        bckInt2, bckBins2 = BinHistograms(Bckdata[1][1],BinRange[1],ChBins[1])
        bckInt4, bckBins4 = BinHistograms(Bckdata[2][1],BinRange[2],ChBins[2])
    
    normalize = norm
    logScale = log
    
    fig, ax = plt.subplots(1,3, figsize = (30,10))

    ax[0].hist(Bins0[:-1],bins =Bins0,density = normalize,weights=Int0/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Ch 0 (LSC)', log = logScale)
    if bck:
        ax[0].hist(bckBins0[:-1],bins =bckBins0,density = normalize,weights=bckInt0/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = 'Background Ch 0 (LSC)', log = logScale)

    ax[1].hist(Bins2[:-1],bins =Bins2,density = normalize,weights=Int2/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = 'Ch 2 (LSC)', log = logScale)
    if bck:
        ax[1].hist(bckBins2[:-1],bins =bckBins2,density = normalize,weights=bckInt2/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = 'Background Ch 0 (LSC)', log = logScale)

    ax[2].hist(Bins4[:-1],bins =Bins4,density = normalize,weights=Int4/(data[0][0][-1] - data[0][0][0]),histtype = 'step',label = f'Ch 4 ({gammaDet})', log = logScale)
    if bck: 
        ax[2].hist(bckBins4[:-1],bins =bckBins4,density = normalize,weights=bckInt4/(Bckdata[0][0][-1] - Bckdata[0][0][0]),histtype = 'step',label = f'Background Ch 0 ({gammaDet})', log = logScale)

    ax[0].set_title('Left PMT')
    ax[1].set_title('Right PMT')
    ax[2].set_title(f'{gammaDet}')

    ax[0].set_xlabel('Integral (ADC Channel)')
    ax[1].set_xlabel('Integral (ADC Channel)')
    ax[2].set_xlabel('Integral (ADC Channel)')

    ax[0].set_ylabel('Counts/bin/ps')
    ax[1].set_ylabel('Counts/bin/ps')
    ax[2].set_ylabel('Counts/bin/ps')
    
    ax[0].legend(loc = 'best')
    ax[1].legend(loc = 'best')
    ax[2].legend(loc = 'best')
    
    
    # plt.show()
    plt.savefig(f'{saveFilePath}/{fileName}_integral_1D_hist.png',bbox_inches='tight')
    # plt.show()
    plt.close()
    
    
    
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
    
    
filepath = '/home/nick/PhD/KDK+/Large_LSC_testing/Position_tests/2025_05_01/'

# GammaDetector = 'LYSO'
# file = 'SDataR_Small_LSC_vessel_LYSO_Cs137_triple_coinc_Vertical_scatter_v4_10CG_LSC.CSV'
# bckfile = '/home/nick/PhD/KDK+/Small_LSC_testing/2025_03_19/SDataR_Small_LSC_vessel_LYSO_bck_triple_coinc_Vertical_scatter_v4_10CG_LSC.CSV'
# # bckfile = f'{filepath}/{file}'

# fileData = ReadInFile(f'{filepath}/{file}', 3)
# filebck = ReadInFile(f'{bckfile}', 3)

GammaDetector = 'NaI'
file = 'SDataR_Large_LSC_vessel_Cs137_triple_coinc_Vertical_scatter_v4_position_3_10fc_corrected.CSV'
bckfile = '/home/nick/PhD/KDK+/Small_LSC_testing/2025_03_19/SDataR_Small_LSC_vessel_NaI_bck_triple_coinc_Vertical_scatter_v4_10CG_LSC.CSV'
# bckfile = f'{filepath}/{file}'

fileData = ReadInFileNaI(f'{filepath}/{file}', 3)
filebck = ReadInFileNaI(f'{bckfile}', 3)

Channels = [0,2,6]

print(f"Ch {Channels[0]}: {len(fileData[0][0])}")
print(f"Ch {Channels[1]}: {len(fileData[1][0])}")
print(f"Ch {Channels[2]}: {len(fileData[2][0])}")

backgrounds = False
dataFile = fileData
bckdataFile = filebck

fileName = file.split('.CSV')[0]
bckfileName = bckfile.split('/')[-1].split('.CSV')[0]

saveFilePath = f"{filepath}/{fileName}/figures"
Path(f"{saveFilePath}").mkdir(parents=True, exist_ok=True)

print(f'Event Rate: {len(dataFile[0][0])/((dataFile[0][0][-1] - dataFile[0][0][0])*1e-12)}')


energyHist1D(data = dataFile, Bckdata = False, ChBins = [100,100,100], BinRange = [[0,4000],[0,4000],[0,4000]], norm = False, log = True,saveFilePath = saveFilePath,fileName = fileName, gammaDet=GammaDetector, bck = backgrounds)
TimeDiff1D(data = dataFile, BinRange = [-300000,300000],saveFilePath = saveFilePath,fileName = fileName)
energyHist2D(dataFile, ChBins = [100,100,100], BinRange = [[0,4000],[0,4000],[0,4000]],saveFilePath = saveFilePath,fileName = fileName)
TimDiff2D(data = dataFile,ChBins = [100,100,100],ChBinRange= [[0,4000],[0,4000],[0,4000]],TimeBinRange= [-400000,400000],saveFilePath = saveFilePath,fileName = fileName)

# energyHist2D(bckdataFile, ChBins = [100,100,100], BinRange = [[0,4000],[0,4000],[0,4000]],saveFilePath = saveFilePath,fileName = bckfileName)
# TimDiff2D(data = bckdataFile,ChBins = [100,100,100],ChBinRange= [[0,4000],[0,4000],[0,4000]],TimeBinRange= [-400000,400000],saveFilePath = saveFilePath,fileName = bckfileName)


TimeCutHist(data=dataFile,ChBins = [100,100,100],ChBinRange = [[0,4000],[0,4000],[0,4000]],norm = False,TimeBinRange= [200000,350000],timeCut = [230000,250000],saveFilePath = saveFilePath,fileName= fileName)
