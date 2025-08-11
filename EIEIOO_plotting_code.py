import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from pathlib import Path

matplotlib.rcParams["font.size"] = 15
matplotlib.rcParams["lines.linewidth"] = 3
matplotlib.rcParams["mathtext.default"] = 'regular'
matplotlib.rcParams['lines.markersize'] = 3

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
    
    
def energyHist1DComp(dataC,dataNoC, ChBins,BinRange,norm,log,saveFilePath,fileName,gammaDet):
    
    Int0, Bins0 = BinHistograms(dataC[0][1],BinRange[0],ChBins[0])
    Int2, Bins2 = BinHistograms(dataC[1][1],BinRange[1],ChBins[1])
    Int4, Bins4 = BinHistograms(dataC[2][1],BinRange[2],ChBins[2])

    NoCInt0, NoCBins0 = BinHistograms(dataNoC[0][1],BinRange[0],ChBins[0])
    NoCInt2, NoCBins2 = BinHistograms(dataNoC[1][1],BinRange[1],ChBins[1])
    NoCInt4, NoCBins4 = BinHistograms(dataNoC[2][1],BinRange[2],ChBins[2])
    
    normalize = norm
    logScale = log
    
    fig, ax = plt.subplots(1,3,figsize = (30,8))
    plt.tight_layout()
    plt.subplots_adjust(left = 0.06,wspace = 0.2)
    ax[0].hist(Bins0[:-1],bins =Bins0,density = normalize,weights=Int0/(dataC[0][0][-1] - dataC[0][0][0]),histtype = 'step',label = 'LSC Left PMT Coincidence', log = logScale)
    ax[0].hist(NoCBins0[:-1],bins =NoCBins0,density = normalize,weights=NoCInt0/(dataNoC[0][0][-1] - dataNoC[0][0][0]),histtype = 'step',label = 'LSC Left PMT No Coincidence', log = logScale)
    
    ax[1].hist(Bins2[:-1],bins =Bins2,density = normalize,weights=Int2/(dataC[0][0][-1] - dataC[0][0][0]),histtype = 'step',label = 'LSC Right PMT Coincidence', log = logScale)
    ax[1].hist(NoCBins2[:-1],bins =NoCBins2,density = normalize,weights=NoCInt2/(dataNoC[0][0][-1] - dataNoC[0][0][0]),histtype = 'step',label = 'LSC Right PMT No Coincidence', log = logScale)
    
    ax[2].hist(Bins4[:-1],bins =Bins4,density = normalize,weights=Int4/(dataC[0][0][-1] - dataC[0][0][0]),histtype = 'step',label = 'NaI Coincidence', log = logScale)
    ax[2].hist(NoCBins4[:-1],bins =NoCBins4,density = normalize,weights=NoCInt4/(dataNoC[0][0][-1] - dataNoC[0][0][0]),histtype = 'step',label = 'NaI No Coincidence', log = logScale)
    
    
    
    ax[0].set_title('Left PMT')
    ax[1].set_title('Right PMT')
    ax[2].set_title(f'{gammaDet}')

    ax[0].set_xlabel('Integral (ADC Channel)')
    ax[1].set_xlabel('Integral (ADC Channel)')
    ax[2].set_xlabel('Integral (ADC Channel)')

    ax[0].set_ylabel('Counts/bin/ps')
    ax[1].set_ylabel('Counts/bin/ps')
    ax[2].set_ylabel('Counts/bin/ps')
    
    ax[0].legend(loc = 'lower center')
    ax[1].legend(loc = 'lower center')
    ax[2].legend(loc = 'lower center')
    
    if norm:
        plt.savefig(f'{saveFilePath}/{fileName}_integral_1D_hist_normalized.png',bbox_inches='tight')
    else:
        plt.savefig(f'{saveFilePath}/{fileName}_integral_1D_hist.png',bbox_inches='tight')
    plt.show()
    
filepath = '/home/nick/PhD/KDK+/Presentation_data/Coincidence_demo_data/'
file = 'SDataR_Large_LSC_vessel_Cs137_triple_coinc_Vertical_scatter_v4_presentation_2.CSV'
fileNoC = 'SDataR_Large_LSC_vessel_Cs137_no_coinc_Vertical_scatter_v4_presentation_2.CSV'

fileData = ReadInFileNaI(f'{filepath}/{file}', 3)
fileDataNoC = ReadInFileNaI(f'{filepath}/{fileNoC}', 3)
dataFile = fileData
dataFileNoC = fileDataNoC

fileName = file.split('.CSV')[0]
saveFilePath = f"{filepath}/{fileName}/figures"
Path(f"{saveFilePath}").mkdir(parents=True, exist_ok=True)

print(f'Event Rate: {len(dataFile[0][0])/((dataFile[0][0][-1] - dataFile[0][0][0])*1e-12)}')
    
energyHist1DComp(dataC = dataFile, dataNoC = dataFileNoC, ChBins = [100,100,100], BinRange = [[0,4000],[0,4000],[0,4000]], norm = False, log = True,saveFilePath = saveFilePath,fileName = fileName,gammaDet="NaI")